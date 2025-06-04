import os
import json
import pickle
import pysam
import numpy as np
from tqdm import tqdm
import random
from random import sample
from collections import Counter
random.seed(0)

mod_tags = {
    'm6A': ('A', 0, 'a'),
    'psi': ('T', 0, 17802),
    '4sU': ('T', 0, 20480)
}


def get_read_mod_prob(in_read, in_mod, in_quantile):
    if in_read.modified_bases is None:
        return 0.0
    pos_prob = in_read.modified_bases.get(mod_tags[in_mod], [(0, 0.0)])
    return np.quantile([prob for pos, prob in pos_prob], in_quantile) / 255.0


def get_mean_mod_prob(in_read, in_mod):
    if in_read.modified_bases is None:
        return 0.0
    pos_prob = in_read.modified_bases.get(mod_tags[in_mod], [(0, 0.0)])
    return np.mean([prob for pos, prob in pos_prob]) / 255.0


def get_read_phred_score(in_read, in_quantile):
    phred_scores = np.array(in_read.query_qualities.tolist()) / 50.0
    return np.quantile(phred_scores, in_quantile)


def get_in_del_ratio(in_read):
    cigartuples = in_read.cigartuples
    op_len = {op: 0 for op in [0, 1, 2]}   # 0: M; 1: I; 2: D
    for tup in cigartuples:
        if tup[0] > 2:
            continue
        else:
            op_len[tup[0]] = op_len[tup[0]] + tup[1]
    total_len = sum(op_len.values())
    ratio_in = op_len[1] / total_len
    ratio_del = op_len[2] / total_len
    return ratio_in + ratio_del


def get_u_ratio(in_read):
    base_counts = Counter(in_read.query_sequence)
    total = np.sum(list(base_counts.values()))
    return base_counts['T'] / total
    # return (base_counts['T'] + base_counts['A']) / (base_counts['C'] + base_counts['G'])


def get_mapq(in_read):
    return in_read.mapping_quality


def get_read_length(in_read, read_len_min=100.0, a_max=5000.0):
    # return np.clip(in_read.query_length, 0, a_max) / a_max
    return read_len_min / np.clip(in_read.query_length, a_min=read_len_min, a_max=None)


def get_mod_occupancy(in_read, in_mod, in_thresh_mod=0.5, min_locs=1):
    if in_read.modified_bases is None:
        mod_mean_occupancy = 0.0
    else:
        read_mod_probs = np.array([this_tup[1] for this_tup in in_read.modified_bases.get(mod_tags[in_mod], [])]) / 255.0
        if len(read_mod_probs) >= min_locs:
            mod_mean_occupancy = np.mean(read_mod_probs >= in_thresh_mod)
        else:
            mod_mean_occupancy = 0.0
    return mod_mean_occupancy


def get_norm_feature_mat(in_feature_mat):
    # out_feature_mat = in_feature_mat - np.mean(in_feature_mat, axis=0)
    # out_feature_mat = in_feature_mat / np.std(in_feature_mat, axis=0)
    out_feature_mat = in_feature_mat / np.max(in_feature_mat, axis=0)
    return out_feature_mat


def get_feature_extractor(in_cfg):
    include_functions = []
    if in_cfg.get('quantile_psi', -1) > 0:
        include_functions.append(lambda x: get_read_mod_prob(x, in_mod='psi', in_quantile=in_cfg['quantile_psi']))
    if in_cfg.get('quantile_m6A', -1) > 0:
        include_functions.append(lambda x: get_read_mod_prob(x, in_mod='m6A', in_quantile=in_cfg['quantile_m6A']))
    if in_cfg.get('quantile_phred', -1) > 0:
        include_functions.append(lambda x: get_read_phred_score(x, in_quantile=in_cfg['quantile_phred']))
    if in_cfg.get('use_in_del_ratio', -1) > 0:
        include_functions.append(lambda x: get_in_del_ratio(x))
    if in_cfg.get('use_u_ratio', -1) > 0:
        include_functions.append(lambda x: get_u_ratio(x))
    if in_cfg.get('use_mapq', -1) > 0:
        include_functions.append(lambda x: get_mapq(x))
    if in_cfg.get('read_len_min', -1) > 0:
        include_functions.append(lambda x: get_read_length(x, read_len_min=in_cfg['read_len_min']))
    if in_cfg.get('thresh_mod', -1) > 0:
        include_functions.append(lambda x: get_mod_occupancy(x, in_mod='psi', in_thresh_mod=in_cfg['thresh_mod']))
        include_functions.append(lambda x: get_mod_occupancy(x, in_mod='m6A', in_thresh_mod=in_cfg['thresh_mod']))
    if in_cfg.get('use_4sU_tag', -1) > 0:
        include_functions.append(lambda x: get_mean_mod_prob(x, in_mod='4sU'))

    def _feature_extractor(in_read):
        out_feature = []
        for this_fcn in include_functions:
            out_feature.append(this_fcn(in_read))
        return out_feature

    return _feature_extractor


def get_features_from_reads(in_feature_extractor, in_bam_file, return_read_name=False):
    print(f'Collecting features from {in_bam_file}')
    read_features = []
    read_names = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            read_features.append(in_feature_extractor(read))
            read_names.append(read.query_name)
    feature_mat = np.vstack(read_features)
    if return_read_name:
        return feature_mat, read_names
    else:
        return feature_mat


def get_feature_and_label(in_bam_pos, in_bam_neg, in_cfg, num_samples=100000):
    feature_extractor = get_feature_extractor(in_cfg)
    feature_pos = get_features_from_reads(feature_extractor, in_bam_pos)
    feature_neg = get_features_from_reads(feature_extractor, in_bam_neg)
    actual_num_samples = min(min(feature_neg.shape[0], feature_pos.shape[0]), num_samples)
    # if actual_num_samples < num_samples:
    #     print(f'Actual num. samples used {actual_num_samples}')

    out_X = np.vstack([
        feature_neg[sample(range(feature_neg.shape[0]), actual_num_samples)],
        feature_pos[sample(range(feature_pos.shape[0]), actual_num_samples)]
    ])
    out_y = np.concatenate([np.zeros(actual_num_samples), np.ones(actual_num_samples)])

    if in_cfg.get('normalize', False):
        out_X = get_norm_feature_mat(out_X)

    return out_X, out_y


def load_model(in_args):
    with open(os.path.join(in_args.model_dir, 'config.json'), 'r') as h_cfg:
        out_cfg = json.load(h_cfg)
    with open(os.path.join(in_args.model_dir, 'clf.pkl'), 'rb') as h_clf:
        out_clf = pickle.load(h_clf)
    print(f"{out_cfg['name']} loaded")
    return out_clf, out_cfg
