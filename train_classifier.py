import os
import pysam
import numpy as np
from tqdm import tqdm
from random import sample
from sklearn.linear_model import LogisticRegression
import pickle
import json


def get_read_mod_prob(in_read, in_quantile, mod_key=('T', 1, 17802)):
    if in_read.modified_bases is None:
        return 0.0
    pos_prob = in_read.modified_bases.get(mod_key, [(0, 0.0)])
    return np.quantile([prob for pos, prob in pos_prob], in_quantile) / 255.0


def get_read_phred_score(in_read, in_quantile):
    phred_scores = np.array(in_read.query_qualities.tolist())
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
    return ratio_in, ratio_del


def get_features_from_reads(in_bam_file, in_quantile_phred, in_quantile_psi):
    mod_prob = []
    phred_score = []
    in_del_ratio = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            mod_prob.append(get_read_mod_prob(read, in_quantile_psi))
            phred_score.append(get_read_phred_score(read, in_quantile_phred))
            in_del_ratio.append(get_in_del_ratio(read))
    feature_mat = np.vstack([mod_prob, np.array(phred_score)/50.0, np.vstack(in_del_ratio).sum(axis=1)]).T
    return feature_mat


def get_feature_and_label(bam_files, num_samples, in_cfg):
    tp_feature = {tp: {} for tp in in_cfg['tps']}
    for tp in in_cfg['tps']:
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


### config ###
cfg = {
    'data_dir': '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU',
    'train_ds': 'chr1',
    'tps': ['0h', '24h'],
    'tp_label': {'0h': 0, '24h': 1},
    'quantile_phred': 0.1,
    'quantile_psi': 0.9,
    'num_train_samples': 100000
}
##############
train_bam_files = {tp: os.path.join(cfg['data_dir'], f"hiPSC-CM_{tp}_4sU_{cfg['train_ds']}.bam") for tp in cfg['tps']}
train_X, train_y = get_feature_and_label(train_bam_files, cfg['num_train_samples'])
clf = LogisticRegression(random_state=0, verbose=True).fit(train_X, train_y)
train_acc = clf.score(train_X, train_y)
print(f"Train on {cfg['train_ds']}, {len(train_y)} reads, accuracy {train_acc:.3f}")

### test ###
test_ds = 'chr2'
num_test_samples = 100000
test_bam_files = {tp: os.path.join(cfg['data_dir'], f'hiPSC-CM_{tp}_4sU_{test_ds}.bam') for tp in cfg['tps']}
test_X, test_y = get_feature_and_label(test_bam_files, num_test_samples, cfg)
test_acc = clf.score(test_X, test_y)
print(f'Test on {test_ds}, {len(test_y)} reads, accuracy {test_acc:.3f}')

### save ###
clf_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/classifier'
clf_name = 'trial0'
os.makedirs(os.path.join(clf_dir, clf_name), exist_ok=True)
with open(os.path.join(clf_dir, clf_name, 'clf.pkl'), 'wb') as out_pkl:
    pickle.dump(clf, out_pkl)
with open(os.path.join(clf_dir, clf_name, 'config.json'), 'w') as out_cfg:
    json.dump(cfg, out_cfg)
