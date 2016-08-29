#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source activate crispr_scr

python $DIR/../../yltkv2.py process_fastq \
       -f $DIR/process_fastq_test.fastq \
       -b 1:5,12:17
       -s 18: \
       -i 6:11 \
       --outfastq AAGCTA,CGTGAT,GCGGAC:$DIR/sample1.fastq;ACATCG,CTGATC,CACTGT:$DIR/sample2.fastq

source deactivate

# AAGCTA    119
# CGTGAT    87
# GCGGAC    67
# ACATCG    67
# CTGATC    56
# CACTGT    52
# Others    52
