import os
import numpy as np
from helper_functions import load_model, get_feature_extractor, get_features_from_reads
from argparse import ArgumentParser


def save_inference_results(in_read_names, in_results, in_args):
    out_dir = os.path.dirname(in_args.out_file)
    if len(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(in_args.out_file, 'w') as f_out:
        for this_read_name, this_result in zip(in_read_names, in_results):
            f_out.write(f'{this_read_name}\t{this_result}\n')
    print(f'Inference results on {len(in_results)} reads written to {in_args.out_file}')


def main():
    parser = ArgumentParser()
    parser.add_argument('--model_dir', type=str, required=True,
                        help='directory containing pickled classifier and config')
    parser.add_argument('--bam', type=str, required=True,
                        help='bam file on which to perform inference')
    parser.add_argument('--out_file', required=True,
                        help='output directory for inference results')
    args = parser.parse_args()

    clf, cfg = load_model(args)
    feature_extractor = get_feature_extractor(cfg)
    test_X, test_read_names = get_features_from_reads(feature_extractor, args.bam, return_read_name=True)
    test_results = clf.predict(test_X)
    test_results = np.int64(test_results)

    save_inference_results(test_read_names, test_results, args)


if __name__ == '__main__':
    main()
