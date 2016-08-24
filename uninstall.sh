#!/bin/bash
# Stop on error
set -e

conda remove --name crispr_scr --all --yes
conda clean -a -y

echo === The conda environment for crispr_scr is successfully removed. ===
