#!/bin/bash
# Stop on error
set -e

CUR_DIR=$(pwd)

conda create -n crispr_scr --file requirements.txt -y -c defaults -c bioconda -c r

############ install additional packages

source activate crispr_scr

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
mkdir -p $CONDA_EXTRA

source deactivate

echo === Installing crispr_scr successfully done. ===
