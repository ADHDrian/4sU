# A simple classifier for 4sU-labelled reads
Package requirements (version numbers not critical):
```
matplotlib==3.10.3
numpy==1.23.5
pandas==2.2.3
pysam==0.22.0
scikit_learn==1.4.2
tqdm==4.64.1
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
Finished
```

The trained model and its config are saved to `${out_dir}`.

## Validate
Validate a classifier on another set of reads with known labels:
```
python3 validate_classifier.py \
--model_dir ./assets/train/svm2 \
--bam_positive /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr3.thresh.bam \
--bam_negative /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr3.thresh.bam
```
The output should look like:
```
Validate on 32290 reads, accuracy 0.875
Finished
```

## Inference
Run inference on reads without known labels:
```
python3 inference.py \
--model_dir ./assets/train/svm2 \
--bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr2.thresh.bam \
--out_file ${out_file}
```
`${out_file}` would be a tsv file with read names and 0/1 (unlabelled/labelled) that looks like `assets/inference/svm_chr2_24h.tsv`.

## Data preparation
The input data that are used for training and testing above are not all the possible reads, but those aligned to the fast-decaying transcripts / genes (in top 50% quantile). A script is provided to filter the raw bam files.
```
python3 filter_reads_by_transcript_decay_rate.py \
--in_bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr1.bam \
--out_bam /prj/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr1.thresh0.5.bam \
--decay_rate_file ./assets/data/gene-estimates-annotated-pulseRTc-0-1-2-4-6-8-16.csv \
--gene_bed_file ./assets/data/gene.ensembl_havana.GRCh38.102.bed
```
Use `--decay_rate_quantile` to adjust the decay rate threshold (default 0.5).
