#!/bin/bash

source activate crispr_scr

python libCheck.py -i sample1.fastq sample2.fastq \
       -t /users/test/CRISPR/Oligos.txt \
       -r 1:9 \
       -g 31:50 \
       -o /users/test/CRISPR/ \
       --picard /users/test/CRISPR/picard.jar

source deactivate
