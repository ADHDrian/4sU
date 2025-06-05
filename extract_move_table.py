import os
import pysam
import numpy as np
from tqdm import tqdm
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt


def get_base_dwell_times(in_bam_file, in_base='T'):
    all_dwell_times = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            stride = read.get_tag('mv')[0]
            mv_tag = np.array(read.get_tag('mv')[1:])
            seq = read.get_forward_sequence()
            base_occupancy = np.cumsum(mv_tag)
            u_locs = np.where(np.array(list(seq)) == in_base)[0]
            dwell_times = []
            for this_u_loc in u_locs:
                dwell_times.append(np.sum(base_occupancy == this_u_loc)*stride)
            all_dwell_times.extend(dwell_times)
    return all_dwell_times


def get_mean_base_dwell_times(in_bam_file, in_base='T', in_max=200.0):
    mean_dwell_times = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            stride = read.get_tag('mv')[0]
            mv_tag = np.array(read.get_tag('mv')[1:])
            seq = read.get_forward_sequence()
            base_occupancy = np.cumsum(mv_tag)
            u_locs = np.where(np.array(list(seq)) == in_base)[0]
            dwell_times = []
            for this_u_loc in u_locs:
                this_dwell_time = np.sum(base_occupancy == this_u_loc)*stride
                if this_dwell_time > 0:
                    dwell_times.append(this_dwell_time)
            # mean_dwell_times.append(np.quantile(dwell_times, q=in_quantile))
            mean_dwell_times.append(np.mean(dwell_times)/in_max)
    return mean_dwell_times


chrom = 'chr1'
tps = ['0h', '24h']
bam_file = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_{}_4sU.{}.thresh.bam'
img_out = '/home/adrian/img_out/4sU'

# tp_dwell_times = {}
# for this_tp in tps:
#     tp_dwell_times[this_tp] = get_base_dwell_times(bam_file.format(this_tp, chrom))

tp_mean_dwell_times = {}
for this_tp in tps:
    tp_mean_dwell_times[this_tp] = get_mean_base_dwell_times(bam_file.format(this_tp, chrom))

bin_max = 200
num_bins = 40
bin_edges = np.linspace(0, bin_max, num_bins+1)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
tp_hist = {}
tp_total_reads = {}
for this_tp in tps:
    this_tp_counts = np.histogram(tp_mean_dwell_times[this_tp], bins=bin_edges)[0]
    num_total_bases = this_tp_counts.sum()
    tp_hist[this_tp] = this_tp_counts / num_total_bases
    tp_total_reads[this_tp] = num_total_bases

tp_colors = {
    '0h': 'b',
    '24h': 'r'
}
plt.figure(figsize=(5, 5))
# plt.subplot(1, 2, 1)
for this_tp in tps:
    # plt.hist(tp_hist[this_tp], label=this_tp)
    # plt.semilogy(bin_centers, tp_hist[this_tp], c=tp_colors[this_tp], label=f'{this_tp} ({tp_total_reads[this_tp]})')
    plt.semilogy(bin_centers, tp_hist[this_tp], c=tp_colors[this_tp], label=f'{this_tp} ({tp_total_reads[this_tp]})')
plt.legend(loc='upper right')
plt.xlabel('Mean base U dwell time per read (signal blocks)')
plt.ylabel('Density')
plt.title(chrom)
# plt.subplot(1, 2, 2)
# for this_tp in tps:
#     plt.semilogy(bin_centers, np.cumsum(tp_hist[this_tp]), c=tp_colors[this_tp], label=this_tp)
# plt.xlim([0, 1000])
# plt.xlabel('Base U dwell time (signal blocks)')
# plt.ylabel('CDF')
plt.savefig(os.path.join(img_out, 'hist_U_dwell_time_per_read.png'), bbox_inches='tight')