import os
import numpy as np

data_dir = '/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/chunks_0h_chr1'

kmer_table_file = os.path.join(data_dir, 'kmer_table.npy')
kmer_table = np.load(kmer_table_file)

sequence_file = os.path.join(data_dir, 'sequence.npy')
sequence = np.load(sequence_file, allow_pickle=True)
