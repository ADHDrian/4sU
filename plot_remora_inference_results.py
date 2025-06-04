import pysam
import numpy as np
import os
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt


def get_mod_prob_per_read(in_bam_file, mod_code = ('T', 0, 20480)):
    mod_prob_per_read = []
    with pysam.AlignmentFile(in_bam_file, 'rb') as bam:
        for this_read in bam.fetch():
            mod_probs = [tup[1] for tup in this_read.modified_bases.get(mod_code, [(0, 0)])]
            mod_prob_per_read.append(np.mean(mod_probs))
    return np.array(mod_prob_per_read) / 255.0


tps = ['0h', '24h']
chrom = 'chr2'
bam_files = {
    this_tp: f'/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/inference/infer_{this_tp}.{chrom}.sorted.bam'
    for this_tp in tps
}

tp_mod_prob_dist = {}
for this_tp, this_bam_file in bam_files.items():
    tp_mod_prob_dist[this_tp] = get_mod_prob_per_read(this_bam_file)


img_out = '/home/adrian/img_out/4sU'
tp_colors = {
    '0h': 'r',
    '24h': 'b'
}

xlim = [0.3, 0.7]

plt.figure(figsize=(4, 4))
for this_tp in tps:
    this_dist = tp_mod_prob_dist[this_tp]
    num_total = len(this_dist)
    num_tagged = np.sum(this_dist > 0)
    plt.hist(
        this_dist[this_dist > 0], bins=100, range=[0, 1],
        alpha=0.5, fc=tp_colors[this_tp], label=f'{this_tp}:\n{num_tagged} tagged\n{num_total} total'
    )
legend = plt.legend(loc='upper left')
plt.xlim(xlim)
plt.xticks([xlim[0], 0.5, xlim[1]], [xlim[0], 0.5, xlim[1]])
plt.xlabel('Mean 4sU probability per read')
plt.ylabel('Read counts')
plt.title(chrom)
plt.savefig(os.path.join(img_out, f'hist_4sU_dist_{chrom}.png'), bbox_inches='tight')