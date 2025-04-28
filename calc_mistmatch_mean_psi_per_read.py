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
from collections import Counter


def get_read_mean_mod_prob(in_read, in_quantile, mod_key=('T', 1, 17802)):
    if in_read.modified_bases is None:
        return 0.0
    pos_prob = in_read.modified_bases.get(mod_key, [(0, 0.0)])
    # return np.mean([prob for pos, prob in pos_prob]) / 255.0
    return np.quantile([prob for pos, prob in pos_prob], in_quantile) / 255.0


# def get_read_mean_phred_score(in_read, base_span):
#     ref_seq = in_read.get_reference_sequence()
#     dict_ref_query = {ref_loc: read_loc for read_loc, ref_loc in in_read.get_aligned_pairs()}
#     ref_seq_locs = in_read.reference_start + np.where(np.array(list(ref_seq)) == 'T')[0]
#     query_seq_locs = [dict_ref_query[this_ref_loc] for this_ref_loc in ref_seq_locs]
#     sel_query_seq_locs = list(set([x for this_loc in query_seq_locs if this_loc is not None
#                                    for x in range(this_loc-2, this_loc+1)]))
#
#     # query_seq = in_read.query_sequence
#     # print(in_read.flag)
#     # print(Counter([query_seq[l] for l in query_seq_locs if l is not None]))
#
#     phred_scores = np.array(in_read.query_qualities.tolist())
#     sel_phred_scores = phred_scores[sel_query_seq_locs]
#
#     return np.mean(sel_phred_scores)

def get_read_mean_phred_score(in_read, in_quantile):
    # ref_seq = in_read.get_reference_sequence()
    # dict_ref_query = {ref_loc: read_loc for read_loc, ref_loc in in_read.get_aligned_pairs()}
    # ref_seq_locs = in_read.reference_start + np.where(np.array(list(ref_seq)) == 'T')[0]
    # query_seq_locs = [dict_ref_query[this_ref_loc] for this_ref_loc in ref_seq_locs]
    # sel_query_seq_locs = list(set([x for this_loc in query_seq_locs if this_loc is not None
    #                                for x in range(this_loc-2, this_loc+1)]))

    # query_seq = in_read.query_sequence
    # print(in_read.flag)
    # print(Counter([query_seq[l] for l in query_seq_locs if l is not None]))

    phred_scores = np.array(in_read.query_qualities.tolist())
    # sel_phred_scores = phred_scores[sel_query_seq_locs]

    return np.quantile(phred_scores, in_quantile)


def get_read_num_bases_matched(in_read):
    ref_seq = in_read.get_reference_sequence()
    query_seq = in_read.query_sequence

    dict_ref_query = {ref_loc: read_loc for read_loc, ref_loc in in_read.get_aligned_pairs()}
    ref_seq_locs = in_read.reference_start + np.where(np.array(list(ref_seq)) == 'T')[0]
    query_seq_locs = [dict_ref_query[this_ref_loc] for this_ref_loc in ref_seq_locs]
    query_matched_bases = [query_seq[this_loc] for this_loc in query_seq_locs if this_loc is not None]
    num_total = len(query_matched_bases)
    num_matched = sum([this_base == 'T' for this_base in query_matched_bases])
    # num_total = len(query_seq_locs)
    # num_matched = len([this_loc for this_loc in query_seq_locs if this_loc is not None])

    return [num_total, num_matched]


def get_in_del_ratio(in_read):
    cigartuples = in_read.cigartuples
    op_len = {op: 0 for op in [0, 1, 2]}   # 0: M; 1: I; 2: D
    for tup in cigartuples:
        if tup[0] > 2:
            continue
        else:
            op_len[tup[0]] = op_len[tup[0]] + tup[1]
    total_len = sum(op_len.values())
    ratio_in = op_len[1] / total_len
    ratio_del = op_len[2] / total_len
    return ratio_in, ratio_del


def get_stats_from_reads(in_bam_file, in_quantile_phred, in_quantile_psi):
    num_total_matched = []
    mean_mod_prob = []
    mean_phred_score = []
    in_del_ratio = []
    mapping_score = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for read in tqdm(bam.fetch()):
            # num_total_matched.append(get_num_bases_matched(read))
            mean_mod_prob.append(get_read_mean_mod_prob(read, in_quantile_psi))
            mean_phred_score.append(get_read_mean_phred_score(read, in_quantile_phred))
            in_del_ratio.append(get_in_del_ratio(read))
            mapping_score.append(read.mapping_quality)

    if len(num_total_matched):
        arr_total_matched = np.vstack(num_total_matched)
        frac_matched = arr_total_matched[:, 1] / arr_total_matched[:, 0]
    else:
        arr_total_matched = []
        frac_matched = []

    in_ratio, del_ratio = np.vstack(in_del_ratio).T
    return arr_total_matched, frac_matched, np.array(mean_mod_prob), np.array(mean_phred_score), np.array(mapping_score), np.array(in_ratio), np.array(del_ratio)


def get_norm_hist(in_data, hist_range=[0, 1], num_bins=100):
    counts, bin_edges = np.histogram(in_data, range=hist_range, bins=num_bins)
    pdf = counts / np.sum(counts)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    return bin_edges, bin_centers, counts, pdf


def plot_scatter(in_res, var_x, var_y, lim_x, lim_y, label_x, label_y, plot_title, out_filename):
    xticks = np.linspace(*lim_x, 5)
    yticks = np.linspace(*lim_y, 5)
    plt.figure(figsize=(5*cm, 5*cm))
    for tp in in_res.keys():
        plt.plot(in_res[tp][var_x], in_res[tp][var_y], '.', c=tp_colors[tp], markersize=1, markeredgewidth=0, label=tp)
    plt.xlim(lim_x)
    plt.ylim(lim_y)
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.title(plot_title)
    _lgnd = plt.legend()
    _lgnd.legend_handles[0].set_markersize(3)
    _lgnd.legend_handles[1].set_markersize(3)
    plt.savefig(out_filename, **fig_kwargs)


def plot_2d_hist(tp_res, var_x, var_y, lim_x, lim_y, label_x, label_y, plot_title, out_filename, sep_slope, sep_portions):
    xticks = np.linspace(*lim_x, 5)
    yticks = np.linspace(*lim_y, 5)

    num_bins = 50
    mat_z, bins_x, bins_y = np.histogram2d(tp_res[var_x], tp_res[var_y],
                                           range=[lim_x, lim_y], bins=num_bins)

    separator_y = bins_x * sep_slope

    plt.figure(figsize=(5*cm, 5*cm))
    plt.imshow(np.log10(mat_z.T+1), vmax=3, extent=lim_x+lim_y, aspect='auto', origin='lower')

    sep_perc = np.round(sep_portions / np.sum(sep_portions) * 100.0, 1)

    plt.plot(bins_x, separator_y, 'r--')
    plt.text(0.5, separator_y[int(num_bins/2)]+5, f'{sep_portions[0]} ({sep_perc[0]}%)', ha='center', va='bottom', c='r')
    plt.text(0.5, separator_y[int(num_bins/2)]-5, f'{sep_portions[1]} ({sep_perc[1]}%)', ha='center', va='top', c='r')

    plt.xlim(lim_x)
    plt.ylim(lim_y)
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.title(plot_title)
    plt.savefig(out_filename, **fig_kwargs)


tp_colors = {
    '0h': 'blue',
    '4h': 'purple',
    '24h': 'red'
}


data_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU'
ds = 'chr1'
# tps = ['0h', '4h', '24h']
tps = ['0h', '24h']
bam_files = {tp: os.path.join(data_dir, f'hiPSC-CM_{tp}_4sU_{ds}.bam') for tp in tps}

quantile_phred = 0.1
quantile_psi = 0.9

img_out = '/home/adrian/img_out/4sU'
os.makedirs(img_out, exist_ok=True)

results = {tp: {} for tp in tps}
for tp in tps:
    results[tp]['arr_total'], results[tp]['frac'], results[tp]['mod_prob'], results[tp]['phred_score'], results[tp]['map_score'], results[tp]['in_ratio'], results[tp]['del_ratio'] = get_stats_from_reads(bam_files[tp], quantile_phred, quantile_psi)
    results[tp]['total_counts'] = len(results[tp]['frac'])
    results[tp]['bin_edges'], results[tp]['bin_centers'], results[tp]['counts'], results[tp]['pdf'] = get_norm_hist(results[tp]['frac'])

### mapQ ###
# plt.figure(figsize=(5*cm, 5*cm))
# for tp in tps:
#     edges, centers, counts, pdf = get_norm_hist(results[tp]['map_score'], hist_range=[0, 60], num_bins=60)
#     cdf = np.cumsum(pdf)
#     plt.plot(centers, pdf, c=tp_colors[tp], label=tp)
# plt.xlabel('mapQ')
# plt.ylabel('PDF')
# plt.xticks(np.linspace(0, 60, 4))
# plt.xlim([0, 60])
# plt.legend(loc='upper left')
# plt.title(ds)
# plt.savefig(os.path.join(img_out, f'mapQ.{FMT}'), **fig_kwargs)

### phred scores ###
plt.figure(figsize=(5*cm, 5*cm))
for tp in tps:
    edges, centers, counts, pdf = get_norm_hist(results[tp]['phred_score'], hist_range=[0, 50], num_bins=50)
    cdf = np.cumsum(pdf)
    plt.plot(centers, pdf, c=tp_colors[tp], label=tp)
plt.xticks(np.linspace(0, 50, 3))
plt.xlim([0, 50])
plt.xlabel('phred score')
plt.ylabel('PDF')
plt.legend(loc='upper right')
plt.title(f'{ds}, q{int(quantile_phred*100)}')
plt.savefig(os.path.join(img_out, f'phred_score_q{int(quantile_phred*100)}.{FMT}'), **fig_kwargs)

### in-del ratio ###
x_max = 0.05
plt.figure(figsize=(10*cm, 5*cm))
dict_in_del = {'in': 'insertion', 'del': 'deletion'}
for subplot_ind, in_del in enumerate(['in', 'del']):
    plt.subplot(1, 2, subplot_ind+1)
    for tp in tps:
        edges, centers, counts, pdf = get_norm_hist(results[tp][f'{in_del}_ratio'], hist_range=[0, x_max], num_bins=50)
        plt.plot(centers, pdf, c=tp_colors[tp], label=tp)
    plt.xticks(np.linspace(0, x_max, 3))
    plt.xlim([0, x_max])
    plt.xlabel(f'Ratio {dict_in_del[in_del]}')
    plt.ylabel('PDF')
    plt.legend(loc='upper right')
plt.suptitle(f'{ds}')
plt.savefig(os.path.join(img_out, f'ratio_in_del_{ds}.{FMT}'), **fig_kwargs)


x_max = 0.1
plt.figure(figsize=(5*cm, 5*cm))
for tp in tps:
    edges, centers, counts, pdf = get_norm_hist(results[tp][f'in_ratio']+results[tp][f'del_ratio'], hist_range=[0, x_max], num_bins=50)
    plt.plot(centers, pdf, c=tp_colors[tp], label=tp)
plt.xticks(np.linspace(0, x_max, 3))
plt.xlim([0, x_max])
plt.xlabel(f'Ratio indel')
plt.ylabel('PDF')
plt.legend(loc='upper right')
plt.title(f'{ds}')
plt.savefig(os.path.join(img_out, f'ratio_in_del_combined_{ds}.{FMT}'), **fig_kwargs)

### scatter num U vs. % matched ###
# plt.figure(figsize=(5*cm, 5*cm))
# for tp in tps:
#     plt.plot(results[tp]['arr_total'][:, 0], results[tp]['frac'], '.', c=tp_colors[tp], markersize=1, markeredgewidth=0, label=f"{tp} ({results[tp]['total_counts']})")
# plt.xlim([0, 6000])
# plt.ylim([0, 1])
# plt.xlabel('Num. U\'s in ref. seq.')
# plt.ylabel('Fraction U\'s matched per read')
# plt.title(f'{ds}')
# lgnd = plt.legend()
# for ind in range(len(tps)):
#     lgnd.legend_handles[ind].set_markersize(3)
# plt.savefig(os.path.join(img_out, f'scatter_numU_frac_matched_chr1.{FMT}'), **fig_kwargs)

### pdf psi per read ###
plt.figure(figsize=(5*cm, 5*cm))
for tp in tps:
    results[tp]['bin_edges_mod_prob'], results[tp]['bin_centers_mod_prob'], results[tp]['counts_mod_prob'], results[tp]['pdf_mod_prob'] = get_norm_hist(results[tp]['mod_prob'])
    plt.semilogy(results[tp]['bin_centers_mod_prob'], results[tp]['pdf_mod_prob'], c=tp_colors[tp], label=f"{tp}")
plt.xlim([0, 1])
plt.ylim([0, 0.5])
plt.xlabel('Mean $P(\Psi)$ per read')
plt.ylabel('PDF')
plt.title(f'{ds}, q{int(quantile_psi*100)}')
plt.legend()
plt.savefig(os.path.join(img_out, f'hist_mean_mod_prob_q{int(quantile_psi*100)}.{FMT}'), **fig_kwargs)

# plt.figure(figsize=(5*cm, 5*cm))
# for tp in tps:
#     results[tp]['bin_edges_mod_prob'], results[tp]['bin_centers_mod_prob'], results[tp]['counts_mod_prob'], results[tp]['pdf_mod_prob'] = get_norm_hist(results[tp]['mod_prob'])
#     plt.plot(results[tp]['bin_centers_mod_prob'], results[tp]['pdf_mod_prob'], c=tp_colors[tp], label=f"{tp}")
# plt.xlim([0.0, 1])
# plt.ylim([0, 0.05])
# plt.xlabel('Mean $P(\Psi)$ per read')
# plt.ylabel('PDF')
# plt.title(f'{ds}, q{int(quantile_psi*100)}')
# plt.legend()
# plt.savefig(os.path.join(img_out, f'hist_mean_mod_prob_q{int(quantile_psi*100)}_zoom.{FMT}'), **fig_kwargs)

### scatter % matched vs psi ###
# plot_scatter(results,
#              var_x='mod_prob',
#              var_y='frac',
#              lim_x=[0, 1],
#              lim_y=[0, 1],
#              label_x='Mean $P(\Psi)$ per read',
#              label_y='Fraction U\'s matched per read',
#              plot_title=f'{ds}',
#              out_filename=os.path.join(img_out, f'scatter_mean_mod_prob_frac_matched_chr1.{FMT}')
#              )

### phred score vs psi ###
plot_scatter(results,
             var_x='mod_prob',
             var_y='phred_score',
             lim_x=[0, 1],
             lim_y=[0, 50],
             label_x='Mean $P(\Psi)$ per read',
             label_y='Mean phred score per read',
             plot_title=f'{ds}',
             out_filename=os.path.join(img_out, f'scatter_mean_mod_prob_phred_score_chr1.{FMT}')
             )

separator_slope = 12.5 / 0.5


def get_read_portions_by_separtor(tp_res, var_x, var_y, sep_slope):
    sep_val = tp_res[var_y] / tp_res[var_x]
    mask = tp_res[var_x] > 0
    prop_above = np.sum((sep_val >= sep_slope) * mask)
    prop_below = np.sum((sep_val < sep_slope) * mask)
    return prop_above, prop_below


for tp in tps:
    portions = get_read_portions_by_separtor(results[tp], 'mod_prob', 'phred_score', separator_slope)

    plot_2d_hist(results[tp],
                 var_x='mod_prob',
                 var_y='phred_score',
                 lim_x=[0, 1],
                 lim_y=[0, 50],
                 label_x=f'Mean $P(\Psi)$ per read, q{int(quantile_psi*100)}',
                 label_y=f'Mean phred score per read, q{int(quantile_phred*100)}',
                 plot_title=f'{ds}, {tp}',
                 out_filename=os.path.join(img_out, f'hist2d_mean_mod_prob_phred_score_chr1_{tp}.{FMT}'),
                 sep_slope=separator_slope,
                 sep_portions=portions
                 )

### CDF frac matched ###
# thresh_min_num_U = 100
# filtered_results = {tp: {} for tp in tps}
# for tp in tps:
#     filtered_results[tp]['frac'] = results[tp]['frac'][results[tp]['arr_total'][:, 0] >= thresh_min_num_U]
#     filtered_results[tp]['mod_prob'] = results[tp]['mod_prob'][results[tp]['arr_total'][:, 0] >= thresh_min_num_U]
#     filtered_results[tp]['num_reads'] = len(filtered_results[tp]['frac'])
#     filtered_results[tp]['bin_edges'], filtered_results[tp]['bin_centers'], filtered_results[tp]['counts'], filtered_results[tp]['pdf'] = get_norm_hist(filtered_results[tp]['frac'])
#     filtered_results[tp]['cdf'] = np.cumsum(filtered_results[tp]['pdf'])
#
# plot_scatter_mod_prob_vs_frac_matched(filtered_results,
#                                       f'{ds}\nMin. {thresh_min_num_U} U\'s in ref. seq.',
#                                       os.path.join(img_out, f'scatter_mean_mod_prob_frac_matched_chr1_thresh{thresh_min_num_U}.{FMT}')
#                                       )
#
# plt.figure(figsize=(5*cm, 5*cm))
# # plt.hist(frac_0h, range=[0, 1], bins=50, log=True, label='0h')
# # plt.hist(frac_24h, range=[0, 1], bins=50, log=True, label='24h')
# # plt.semilogy(bin_centers_0h, pdf_0h, c='b', label=f'0h filtered ({num_reads_0h})')
# # plt.semilogy(bin_centers_24h, pdf_24h, c='r', label=f'24h filtered ({num_reads_24h})')
# for tp in tps:
#     plt.semilogy(filtered_results[tp]['bin_centers'], filtered_results[tp]['cdf'], c=tp_colors[tp], label=f"{tp} ({filtered_results[tp]['num_reads']})")
# # plt.plot(bin_centers_0h, cdf_0h, c='b', label=f'0h ({num_reads_0h})')
# # plt.plot(bin_centers_24h, cdf_24h, c='r', label=f'24h ({num_reads_24h})')
# plt.xlim([0, 1])
# plt.legend(loc='lower right')
# plt.title(f'{ds}\nMin. {thresh_min_num_U} U\'s in ref. seq.')
# plt.xlabel('Fraction U\'s matched per read')
# plt.ylabel('CDF reads')
# plt.savefig(os.path.join(img_out, f'cdf_frac_matched_chr1_thresh{thresh_min_num_U}.{FMT}'), **fig_kwargs)
