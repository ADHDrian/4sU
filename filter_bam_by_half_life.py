import os
import pandas as pd
import numpy as np
import pysam
from tqdm import tqdm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def filter_bam_by_gene_bed(in_tp, in_chrom, in_df_gene_bed):
    in_bam_file = os.path.join(bam_dir, f'hiPSC-CM_{in_tp}_4sU_chr{in_chrom}.bam')
    out_bam_file = os.path.join(bam_dir, f'hiPSC-CM_{in_tp}_4sU_chr{in_chrom}.thresh.bam')

    print(f'Parsing {len(in_df_gene_bed)} genes')

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


def get_thresholded_df_gene_bed(in_rate_quantile):
    df_rate = pd.read_csv(half_life_file)

    # xlim = [0, 2]
    # plt.figure(figsize=(5, 5))
    # plt.hist(df_rate['rate'], bins=50, range=xlim)
    # plt.axvline(x=thresh_rate, c='r')
    # plt.xlim(xlim)
    # plt.xlabel('Decay rate')
    # plt.ylabel('Gene count')
    # plt.savefig(os.path.join(img_out, 'hist_rate.png'), bbox_inches='tight')

    thresh_rate = np.quantile(df_rate['rate'], in_rate_quantile)
    df_rate_thresh = df_rate[df_rate['rate'] >= thresh_rate]
    df_gene_bed = pd.read_csv(gene_bed_file, sep='\t')
    out_df_gene_bed = df_gene_bed[
        df_gene_bed['gene_id'].isin(df_rate_thresh['ensembl_gene_id'])
        * df_gene_bed['chrom'] == this_chrom
        ]
    return out_df_gene_bed


bam_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/'
half_life_file = '/home/adrian/Data/TRR319_RMaP_B01/gene-estimates-annotated-pulseRTc-0-1-2-4-6-8-16.csv'
gene_bed_file = '~/Data/genomes/homo_sapiens/GRCh38_102/gene.ensembl_havana.GRCh38.102.bed'
img_out = '/home/adrian/img_out/4sU'

this_chrom = '4'
rate_quantile = 0.5

df_gene_bed_thresh = get_thresholded_df_gene_bed(rate_quantile)

for this_tp in ['0h', '24h']:
    filter_bam_by_gene_bed(this_tp, this_chrom, df_gene_bed_thresh)