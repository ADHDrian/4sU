import pysam
import numpy as np
from tqdm import tqdm
from random import sample
from collections import Counter


def get_read_mod_prob(in_read, in_quantile, mod_key=('T', 1, 17802)):
    if in_read.modified_bases is None:
        return 0.0
    pos_prob = in_read.modified_bases.get(mod_key, [(0, 0.0)])
    return np.quantile([prob for pos, prob in pos_prob], in_quantile) / 255.0


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
    # total = np.sum(list(base_counts.values()))
    return (base_counts['T'] + base_counts['A']) / (base_counts['C'] + base_counts['G'])


def get_read_length(in_read, a_max=5000):
    # return np.clip(in_read.query_length, 0, a_max) / a_max
    return 1.0 / in_read.query_length


def get_features_from_reads(in_bam_file, in_quantile_phred, in_quantile_psi):
    read_features = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            read_features.append([
                get_read_mod_prob(read, in_quantile_psi),
                get_read_phred_score(read, in_quantile_phred),
                get_in_del_ratio(read)
            ])
    feature_mat = np.vstack(read_features)
    return feature_mat


def get_feature_and_label(bam_files, num_samples, in_cfg):
    tp_feature = {tp: {} for tp in in_cfg['tps']}
    for tp in in_cfg['tps']:
        print(f'Collecting features from {tp}:')
        tp_feature[tp] = get_features_from_reads(bam_files[tp], in_cfg['quantile_phred'], in_cfg['quantile_psi'])

    actual_num_samples = min(min(tp_feature['0h'].shape[0], tp_feature['24h'].shape[0]), num_samples)
    if actual_num_samples < num_samples:
        print(f'Actual num. samples used {actual_num_samples}')

    out_X = np.vstack([
        tp_feature['0h'][sample(range(tp_feature['0h'].shape[0]), actual_num_samples)],
        tp_feature['24h'][sample(range(tp_feature['24h'].shape[0]), actual_num_samples)]
    ])
    out_y = np.concatenate([np.zeros(actual_num_samples), np.ones(actual_num_samples)])

    return out_X, out_y