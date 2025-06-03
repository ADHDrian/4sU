# A simple classifier for 4sU-labelled reads
Tested with python3.10. To set up a virtual environment:
```
module load python3/3.10.13_deb12
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

## Train
```
python3 train_classifier.py \
--bam_positive /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr1.thresh.bam \
--bam_negative /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr1.thresh.bam \
--config ./assets/configs/svm2.json \
--out_dir ${out_dir}
```

- `--bam_positive` contains reads that are labelled, `--bam_negative` unlabelled.
- `--config` defines the hyperparameters of the classifier.

When the training is finished (should take less than 10 minutes), it will displaying the following:
```
svm2
Trained on 60000 reads, accuracy 0.794
```

The trained model and its config are saved to `${out_dir}` and should be similar to `assets/train/svm2`.

## Validate
Validate a classifier on another set of reads with known labels:
```
python3 validate_classifier.py \
--model_dir ${out_dir}/svm2 \
--bam_positive /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr3.thresh.bam \
--bam_negative /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr3.thresh.bam
```
The output should look like:
```
Validate on 32290 reads, accuracy 0.875
```

## Inference
Run inference on reads without known labels:
```
python3 inference.py \
--model_dir ${out_dir}/svm2 \
--bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr2.thresh.bam \
--out_file ${out_file}
```
`${out_file}` would be a tsv file with read names and 0/1 (unlabelled/labelled) that looks like `assets/inference/svm_chr2_24h.tsv`.

## Data preparation
The input data that are used for training and testing above are not all the reads, but those aligned to the fast-decaying transcripts / genes (in top 50% quantile). A script is provided to filter the raw bam files.
```
python3 filter_reads_by_transcript_decay_rate.py \
--in_bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr1.bam \
--out_bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr1.thresh.bam \
--decay_rate_file ./assets/data/gene-estimates-annotated-pulseRTc-0-1-2-4-6-8-16.csv \
--gene_bed_file ./assets/data/gene.ensembl_havana.GRCh38.102.bed
```
Use `--decay_rate_quantile` to adjust the decay rate quantile threshold (default 0.5, ie, q50).

## Features
The default configuration uses the following features from each read:
```
quantile_phred	:	0.1   # q10 value in phred score distribution
quantile_psi	:	0.99   # q99 value in P(psi) distribution
use_in_del_ratio	:	1   # (I + D) / (M + I + D)
use_u_ratio	:	1   # num(U) / num(all bases)
use_mapq	:	1   # read mapping quality
read_len_min	:	50   # read_len_min / np.clip(in_read.query_length, a_min=read_len_min)
thresh_mod	:	0.5   # probability threshold for counting an A as m6A / U as psi
```
To turn off a feature, set value to -1 and retrain.

# Note on using remora, kmer and move tables
Besides using read-level features, we also experimented with using ONT remora (https://github.com/nanoporetech/remora) together with their RNA004 kmer models (https://github.com/nanoporetech/kmer_models) to predict position-based 4sU tags directly from the raw signal and move table. The following is a sketch of the steps taken:
## Data preparation
Recall the pod5 reads with `--emit-moves` turned on in dorado:
```
module load dorado/0.9.1

pod5=/prj/TRR319_RMaP_B01/Isabel/20240702_hiPSC-CM_4sU_RNA004/hiPSC-CM_0h_4sU_RTA/20240702_1607_1E_PAS90949_18e6fa36/pod5
bam=/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_0h_4sU_RTA.bam
model=/home/achan/dorado_models/rna004_130bps_sup@v5.1.0

sbatch -p gpu -c 4 --mem 100GB \
--wrap="dorado basecaller ${model} ${pod5} --recursive --modified-bases m5C inosine_m6A pseU -x cuda:0 --estimate-poly-a --emit-moves > ${bam}"
```
Align the reads, threshold them by decay rate as before, then use remora to chunkify the raw signals. The `--reverse-signal` flag is essential for dRNA, as the signal starts from the 3' end, unlike the DNA signal that starts from 5'.
```
pod5="/prj/TRR319_RMaP_B01/Isabel/20240702_hiPSC-CM_4sU_RNA004/hiPSC-CM_24h_4sU_RTA/20240702_1607_1F_PAS90977_9ece4d04/pod5"
bam="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_24h_4sU.chr1.thresh.bam"
levels="/home/achan/git/kmer_models/rna004/9mer_levels_v1.txt"
out_dir="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/chunks_24h_chr1"

srun -c 40 --mem 150GB \
remora \
  dataset prepare \
  ${pod5} \
  ${bam} \
  --output-path ${out_dir} \
  --refine-kmer-level-table ${levels} \
  --refine-rough-rescale \
  --reverse-signal \
  --motif T 0 \
  --mod-base 20480 4sU
  

pod5="/prj/TRR319_RMaP_B01/Isabel/20240702_hiPSC-CM_4sU_RNA004/hiPSC-CM_0h_4sU_RTA/20240702_1607_1E_PAS90949_18e6fa36/pod5"
bam="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_0h_4sU.chr1.thresh.bam"
levels="/home/achan/git/kmer_models/rna004/9mer_levels_v1.txt"
out_dir="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/chunks_0h_chr1"

srun -c 40 --mem 150GB \
remora \
  dataset prepare \
  ${pod5} \
  ${bam} \
  --output-path ${out_dir} \
  --refine-kmer-level-table ${levels} \
  --refine-rough-rescale \
  --reverse-signal \
  --motif T 0 \
  --mod-base-control
```
Create a json config for the datasets:
```
json="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/train_dataset.json"
ctrl="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/chunks_0h_chr1"
mod="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/chunks_24h_chr1"
log="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/train_dataset.log"

srun -c 40 --mem 150GB \
remora \
  dataset make_config \
  ${json} \
  ${ctrl} \
  ${mod} \
  --dataset-weights 1 1 \
  --log-filename ${log}
```
## Training
Once the datasets are prepared, the rest just follows the remora documentation:
```
model="/home/achan/git/remora/models/ConvLSTM_w_ref.py"
json="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/train_dataset.json"
results="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/train_results"

srun -p gpu --gres=gpu:hopper:1 -c 10 --mem 150GB \
remora \
  model train \
  ${json} \
  --model ${model} \
  --device 0 \
  --chunk-context 50 50 \
  --output-path ${results}
```
## Inference
To run a trained model on other datasets:
```
model="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/train_results/model_best.pt"

tp="0h"
chrom="chr3"

pod5="/prj/TRR319_RMaP_B01/Isabel/20240702_hiPSC-CM_4sU_RNA004/hiPSC-CM_${tp}_4sU_RTA/*/pod5"
in_bam="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/hiPSC-CM_${tp}_4sU.${chrom}.thresh.bam"
out_bam="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/inference/infer_${tp}.${chrom}.bam"
log="/prj/TRR319_RMaP_B01/Adrian/4sU/move_table/remora/inference/infer_${tp}.${chrom}.log"

srun -p gpu --gres=gpu:hopper:1 -c 10 --mem 150GB \
remora \
  infer from_pod5_and_bam \
  ${pod5} \
  ${in_bam} \
  --model ${model} \
  --reference-anchored \
  --out-bam ${out_bam} \
  --log-filename ${log} \
  --device 0
```
