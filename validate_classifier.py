from helper_functions import load_model, get_feature_and_label
from argparse import ArgumentParser


def main():
    parser = ArgumentParser()
    parser.add_argument('--model_dir', type=str, required=True,
                        help='directory containing pickled classifier and config')
    parser.add_argument('--bam_positive', type=str, required=True,
                        help='bam file containing positive samples')
    parser.add_argument('--bam_negative', type=str, required=True,
                        help='bam file containing negative samples')
    args = parser.parse_args()

    clf, cfg = load_model(args)
    validate_X, validate_y = get_feature_and_label(args.bam_positive, args.bam_negative, cfg)
    validate_acc = clf.score(validate_X, validate_y)
    print(f'Validate on {len(validate_y)} reads, accuracy {validate_acc:.3f}')


if __name__ == '__main__':
    main()
    print('Finished')
