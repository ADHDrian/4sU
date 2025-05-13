import os
import json
import pickle
from helper_functions import get_feature_and_label


# load #
clf_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/classifier'
clf_name = 'trial11'
with open(os.path.join(clf_dir, clf_name, 'config.json'), 'r') as h_cfg:
    cfg = json.load(h_cfg)
with open(os.path.join(clf_dir, clf_name, 'clf.pkl'), 'rb') as h_clf:
    clf = pickle.load(h_clf)

# test #
test_ds = 'chr2'
num_test_samples = 100000
test_bam_files = {tp: os.path.join(cfg['data_dir'], f'hiPSC-CM_{tp}_4sU_{test_ds}.bam') for tp in cfg['tps']}
test_X, test_y = get_feature_and_label(test_bam_files, num_test_samples, cfg)
test_acc = clf.score(test_X, test_y)
print(f'{clf_name}')
print(f'Test on {test_ds}, {len(test_y)} reads, accuracy {test_acc:.3f}')
