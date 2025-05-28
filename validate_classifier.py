import os
import json
import pickle
from helper_functions import get_feature_and_label
import numpy as np

# load #
clf_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/classifier'
clf_name = 'svm2_thresh'
num_test_samples = 100000

with open(os.path.join(clf_dir, clf_name, 'config.json'), 'r') as h_cfg:
    cfg = json.load(h_cfg)
with open(os.path.join(clf_dir, clf_name, 'clf.pkl'), 'rb') as h_clf:
    clf = pickle.load(h_clf)
print(f'{clf_name}')

# test #
chr_score = {}
for this_chr in [str(this_chr) for this_chr in range(2, 23)] + ['X']:
    test_ds = f'chr{this_chr}'
    test_bam_files = {tp: os.path.join(cfg['data_dir'], f'hiPSC-CM_{tp}_4sU_{test_ds}.thresh.bam') for tp in cfg['tps']}
    # test_bam_files = {tp: os.path.join(cfg['data_dir'], f'hiPSC-CM_{tp}_4sU_{test_ds}.bam') for tp in cfg['tps']}
    test_X, test_y = get_feature_and_label(test_bam_files, num_test_samples, cfg)
    test_acc = clf.score(test_X, test_y)
    chr_score[this_chr] = test_acc
    print(f'Test on {test_ds}, {len(test_y)} reads, accuracy {test_acc:.3f}')

for k, v in chr_score.items():
    print(f'chr{k}: {v:.3f}')
avg_test_acc = np.mean(list(chr_score.values()))
print(f'Mean accuracy: {avg_test_acc:.3f}')
