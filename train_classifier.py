import os
from sklearn.svm import SVC
import pickle
import json
from helper_functions import get_feature_and_label
from argparse import ArgumentParser


def train_model(in_train_X, in_train_y, in_cfg):
    clf = SVC(kernel=in_cfg['svm_kernel'], C=in_cfg['svm_c'], gamma="auto", degree=in_cfg['svm_degree'], verbose=True)
    clf.fit(in_train_X, in_train_y)
    train_acc = clf.score(in_train_X, in_train_y)
    print(in_cfg['name'])
    print(f"Trained on {len(in_train_y)} reads, accuracy {train_acc:.3f}")
    return clf


def save_model(in_clf, in_cfg, in_args):
    os.makedirs(os.path.join(in_args.out_dir, in_cfg['name']), exist_ok=True)
    with open(os.path.join(in_args.out_dir, in_cfg['name'], 'clf.pkl'), 'wb') as out_pkl:
        pickle.dump(in_clf, out_pkl)
    with open(os.path.join(in_args.out_dir, in_cfg['name'], 'config.json'), 'w') as f_out:
        json.dump(in_cfg, f_out)


def main():
    parser = ArgumentParser()
    parser.add_argument('--bam_positive', type=str, required=True,
                        help='bam file containing positive samples')
    parser.add_argument('--bam_negative', type=str, required=True,
                        help='bam file containing negative samples')
    parser.add_argument('--config', type=str, required=True,
                        help='model configuration in json format')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='model output directory')
    args = parser.parse_args()

    with open(args.config, 'r') as f_in:
        cfg = json.load(f_in)

    train_X, train_y = get_feature_and_label(args.bam_positive, args.bam_negative, cfg)
    clf_trained = train_model(train_X, train_y, cfg)
    save_model(clf_trained, cfg, args)


if __name__ == '__main__':
    main()
    print('Finished')
