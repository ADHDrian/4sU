import pysam
from tqdm import tqdm

tp = '24h'
chrom = 'chr3'
donoar_bam_file = f'/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/inference/infer_{tp}.{chrom}.sorted.bam'
acceptor_bam_file = f'/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_{tp}_4sU.{chrom}.thresh.bam'
out_bam_file = f'/home/adrian/Data/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_{tp}_4sU.{chrom}.thresh.4sU.bam'

mod_code = ('T', 0, 20480)

read_tags = {}
with pysam.AlignmentFile(donoar_bam_file, 'rb') as donor_bam:
    for this_read in tqdm(donor_bam.fetch()):
        this_read_mm_tag = this_read.get_tag('MM')
        this_read_ml_tag = this_read.get_tag('ML')
        if (this_read_mm_tag is not None) and (this_read_ml_tag is not None):
            read_tags[this_read.query_name] = (this_read_mm_tag, this_read_ml_tag)

with pysam.AlignmentFile(acceptor_bam_file, 'rb') as acceptor_bam:
    with pysam.AlignmentFile(out_bam_file, 'wb', template=acceptor_bam) as out_bam:
        for this_read in tqdm(acceptor_bam.fetch()):
            if this_read.query_name in read_tags.keys():
                orig_mm = this_read.get_tag('MM')
                orig_ml = this_read.get_tag('ML')
                this_read_tags = read_tags[this_read.query_name]
                new_mm = orig_mm + this_read_tags[0]
                new_ml = orig_ml + this_read_tags[1]
                this_read.set_tag('MM', new_mm)
                this_read.set_tag('ML', new_ml)
            out_bam.write(this_read)
pysam.index(out_bam_file)
