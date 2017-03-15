#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source activate crispr_scr

python $DIR/../../guidescount.py -i $DIR/test1.fastq $DIR/test2.fastq \
       -t $DIR/Table_Guides.txt \
       -r 1:9 \
       -g 31:50 \
       -o $DIR/

source deactivate
