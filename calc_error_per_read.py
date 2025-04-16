import os
import pysam
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
#######################################################################
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Arial'
# FMT = 'svg'
# fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
FMT = 'png'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi)
#######################################################################
import matplotlib.pyplot as plt
from tqdm import tqdm


def get_frac_matched(in_bam_file):
    num_total_matched = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            # print(read.get_aligned_pairs(with_seq=True))

            ref_seq = read.get_reference_sequence()
            query_seq = read.query_sequence

            dict_ref_query = {ref_loc: read_loc for read_loc, ref_loc in read.get_aligned_pairs()}
            ref_seq_locs = read.reference_start + np.where(np.array(list(ref_seq)) == 'T')[0]
            query_seq_locs = [dict_ref_query[this_ref_loc] for this_ref_loc in ref_seq_locs]
            query_matched_bases = [query_seq[this_loc] for this_loc in query_seq_locs if this_loc is not None]

            num_total = len(query_matched_bases)
            num_matched = sum([this_base == 'T' for this_base in query_matched_bases])

            # num_total = len(query_seq_locs)
            # num_matched = len([this_loc for this_loc in query_seq_locs if this_loc is not None])
            num_total_matched.append([num_total, num_matched])

    arr_total_matched = np.vstack(num_total_matched)
    frac_matched = arr_total_matched[:, 1] / arr_total_matched[:, 0]
    return arr_total_matched, frac_matched


def get_norm_hist(in_frac, hist_range=[0, 1], num_bins=100):
    counts, bin_edges = np.histogram(in_frac, range=hist_range, bins=num_bins)
    pdf = counts / np.sum(counts)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return bin_edges, bin_centers, counts, pdf


data_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU'
ds = 'chr1'
bam_file_0h = os.path.join(data_dir, f'hiPSC-CM_0h_4sU_{ds}.bam')
bam_file_24h = os.path.join(data_dir, f'hiPSC-CM_24h_4sU_{ds}.bam')

img_out = '/home/adrian/img_out/4sU'
os.makedirs(img_out, exist_ok=True)

arr_total_0h, frac_0h = get_frac_matched(bam_file_0h)
arr_total_24h, frac_24h = get_frac_matched(bam_file_24h)

total_counts_0h = len(frac_0h)
total_counts_24h = len(frac_24h)

bin_edges_0h, bin_centers_0h, counts_0h, pdf_0h = get_norm_hist(frac_0h)
bin_edges_24h, bin_centers_24h, counts_24h, pdf_24h = get_norm_hist(frac_24h)

plt.figure(figsize=(5*cm, 5*cm))
plt.plot(arr_total_0h[:, 0], frac_0h, '.', markersize=1, markeredgewidth=0, label=f'0h ({total_counts_0h})')
plt.plot(arr_total_24h[:, 0], frac_24h, '.', markersize=1, markeredgewidth=0, label=f'24h ({total_counts_24h})')
plt.xlim([0, 6000])
plt.ylim([0, 1])
plt.xlabel('Num. U\'s in ref. seq.')
plt.ylabel('Fraction U\'s matched per read')
plt.title(f'{ds}')
lgnd = plt.legend()
lgnd.legend_handles[0].set_markersize(3)
lgnd.legend_handles[1].set_markersize(3)
plt.savefig(os.path.join(img_out, f'scatter_numU_frac_matched_chr1.{FMT}'), **fig_kwargs)

thresh_min_num_U = 1000
filtered_frac_0h = frac_0h[arr_total_0h[:, 0] >= thresh_min_num_U]
filtered_frac_24h = frac_24h[arr_total_24h[:, 0] >= thresh_min_num_U]

num_reads_0h = len(filtered_frac_0h)
num_reads_24h = len(filtered_frac_24h)

bin_edges_0h, bin_centers_0h, counts_0h, pdf_0h = get_norm_hist(filtered_frac_0h)
bin_edges_24h, bin_centers_24h, counts_24h, pdf_24h = get_norm_hist(filtered_frac_24h)

cdf_0h = np.cumsum(pdf_0h)
cdf_24h = np.cumsum(pdf_24h)

plt.figure(figsize=(5*cm, 5*cm))
# plt.hist(frac_0h, range=[0, 1], bins=50, log=True, label='0h')
# plt.hist(frac_24h, range=[0, 1], bins=50, log=True, label='24h')
# plt.semilogy(bin_centers_0h, pdf_0h, c='b', label=f'0h filtered ({num_reads_0h})')
# plt.semilogy(bin_centers_24h, pdf_24h, c='r', label=f'24h filtered ({num_reads_24h})')
plt.semilogy(bin_centers_0h, cdf_0h, c='b', label=f'0h ({num_reads_0h})')
plt.semilogy(bin_centers_24h, cdf_24h, c='r', label=f'24h ({num_reads_24h})')
# plt.plot(bin_centers_0h, cdf_0h, c='b', label=f'0h ({num_reads_0h})')
# plt.plot(bin_centers_24h, cdf_24h, c='r', label=f'24h ({num_reads_24h})')
plt.xlim([0, 1])
plt.legend(loc='lower right')
plt.title(f'{ds}\nMin. {thresh_min_num_U} U\'s in ref. seq.')
plt.xlabel('Fraction U\'s matched per read')
plt.ylabel('CDF reads')
plt.savefig(os.path.join(img_out, f'cdf_frac_matched_chr1_thresh{thresh_min_num_U}.{FMT}'), **fig_kwargs)
