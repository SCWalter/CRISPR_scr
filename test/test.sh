#!/bin/bash

source activate crispr_scr

python ../guideCount.py -i test1.fastq test2.fastq \
       -t ./Table_Guides.txt \
       -r 1:9 \
       -g 31:50 \
       -o ./

source deactivate
