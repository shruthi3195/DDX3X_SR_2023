#!/bin/bash

#This runs Kallisto on all paired files in a directory. 

#There are 4 files: all have the same sample number "S1" - in lane 1 L001, in lane 2 L002, fwd read, R1, rev read, R2
for f in *_L001_R1_001.fastq.gz; do
    name=$(basename $f _L001_R1_001.fastq.gz)
    fileL2R1=$(basename $f _L001_R1_001.fastq.gz)_L002_R1_001.fastq.gz
    fileL1R2=$(basename $f _L001_R1_001.fastq.gz)_L001_R2_001.fastq.gz
    fileL2R2=$(basename $f _L001_R1_001.fastq.gz)_L002_R2_001.fastq.gz
    mkdir /lab/page_human_data/shruthir/Azfa/kallisto/$name
    bsub -n 16 "kallisto quant -i index.ksto -t --bootstrap-samples=100 16 --bias --plaintext -o kallisto/$name/ $f $fileL1R2 $fileL2R1 $fileL2R2"

done
