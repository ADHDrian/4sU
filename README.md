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
```
python3 validate_classifier.py \
--model_dir ./assets/train/svm2 \
--bam_positive /prj/TRR319_RMaP_B01/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr3.thresh.bam \
--bam_negative /prj/TRR319_RMaP_B01/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_0h_4sU_chr3.thresh.bam
```

## Inference
```
python3 inference.py \
--model_dir ./assets/train/svm2 \
--bam /prj/TRR319_RMaP_B01/TRR319_RMaP_B01/Adrian/4sU/bam/hiPSC-CM_24h_4sU_chr2.thresh.bam \
--out_file ${out_file}
```
