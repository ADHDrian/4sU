import os
from sklearn.linear_model import LogisticRegression
import pickle
import json
from helper_functions import get_feature_and_label

clf_name = 'trial1'
cfg = {
    'data_dir': '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU',
    'train_ds': 'chr1',
    'tps': ['0h', '24h'],
    'tp_label': {'0h': 0, '24h': 1},
    'quantile_phred': 0.1,
    'quantile_psi': 0.9,
    'num_train_samples': 100000
}

# train #
train_bam_files = {tp: os.path.join(cfg['data_dir'], f"hiPSC-CM_{tp}_4sU_{cfg['train_ds']}.bam") for tp in cfg['tps']}
train_X, train_y = get_feature_and_label(train_bam_files, cfg['num_train_samples'], cfg)
clf = LogisticRegression(random_state=0, verbose=True).fit(train_X, train_y)
train_acc = clf.score(train_X, train_y)
print(f"Train on {cfg['train_ds']}, {len(train_y)} reads, accuracy {train_acc:.3f}")

# save #
clf_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/classifier'
os.makedirs(os.path.join(clf_dir, clf_name), exist_ok=True)
with open(os.path.join(clf_dir, clf_name, 'clf.pkl'), 'wb') as out_pkl:
    pickle.dump(clf, out_pkl)
with open(os.path.join(clf_dir, clf_name, 'config.json'), 'w') as out_cfg:
    json.dump(cfg, out_cfg)
