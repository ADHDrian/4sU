import os
import pandas as pd
import numpy as np
import pysam
from tqdm import tqdm
from argparse import ArgumentParser


def filter_bam_by_gene_bed(in_bam_file, out_bam_file, in_df_gene_bed):
    print(f'Selecting reads in {len(in_df_gene_bed)} genes')

    num_reads_read = 0
    num_reads_written = 0
    with pysam.AlignmentFile(in_bam_file, 'rb') as in_bam:
        with pysam.AlignmentFile(f'{out_bam_file}.temp', 'wb', template=in_bam) as out_bam:
            for _, this_row in tqdm(in_df_gene_bed.iterrows()):
                if this_row['strand'] == '+':
                    flag_required = 0
                elif this_row['strand'] == '-':
                    flag_required = 16
                for this_read in in_bam.fetch(this_row['chrom'], this_row['chromStart'], this_row['chromEnd']):
                    num_reads_read += 1
                    if this_read.flag == flag_required:
                        out_bam.write(this_read)
                        num_reads_written += 1
    print(f'{num_reads_read} reads in, {num_reads_written} out')

    pysam.sort("-o", out_bam_file, f'{out_bam_file}.temp')
    pysam.index(out_bam_file)
    os.remove(f'{out_bam_file}.temp')


def get_thresholded_df_gene_bed(in_decay_rate_file, in_gene_bed_file, in_rate_quantile):
    df_rate = pd.read_csv(in_decay_rate_file)
    thresh_rate = np.quantile(df_rate['rate'], in_rate_quantile)
    df_rate_thresh = df_rate[df_rate['rate'] >= thresh_rate]
    df_gene_bed = pd.read_csv(in_gene_bed_file, sep='\t')
    out_df_gene_bed = df_gene_bed[
        df_gene_bed['gene_id'].isin(df_rate_thresh['ensembl_gene_id'])
        # * df_gene_bed['chrom'] == this_chrom
        ]
    return out_df_gene_bed


def main():
    parser = ArgumentParser()
    parser.add_argument('--in_bam', type=str, required=True,
                        help='input bam file containing original reads')
    parser.add_argument('--out_bam', type=str, required=True,
                        help='output bam file containing reads filtered by transcript decay rate')
    parser.add_argument('--decay_rate_file', type=str, required=True,
                        help='csv file listing gene / transcript decay rate')
    parser.add_argument('--gene_bed_file', type=str, required=True,
                        help='bed file for gene annotation')
    parser.add_argument('--decay_rate_quantile', type=float, default=0.5,
                        help='quantile of decay rate above which reads are kept')
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out_bam)
    if len(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    df_gene_bed_thresh = get_thresholded_df_gene_bed(args.decay_rate_file, args.gene_bed_file, args.decay_rate_quantile)
    filter_bam_by_gene_bed(args.in_bam, args.out_bam, df_gene_bed_thresh)


if __name__ == '__main__':
    main()
