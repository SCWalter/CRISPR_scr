#!/bin/bash
# Stop on error
set -e

conda create -n crispr_scr --file requirements.txt -y -c defaults -c bioconda

############ install additional packages

source activate crispr_scr

CONDA_ENV=$(dirname $(dirname $(which activate)))
CONDA_ACTIVATE_D="$CONDA_ENV/etc/conda/activate.d"
CONDA_DEACTIVATE_D="$CONDA_ENV/etc/conda/deactivate.d"
mkdir -p $CONDA_ACTIVATE_D $CONDA_DEACTIVATE_D

# Picard
touch $CONDA_ACTIVATE_D"/env_vars.sh"
PICARD_PATH=$CONDA_ENV"/share/picard-2.5.0-1/picard.jar"
printf '#!/bin/sh\n\nexport PICARD='"'"$PICARD_PATH"'"'\n' > $CONDA_ACTIVATE_D"/env_vars.sh"
touch $CONDA_DEACTIVATE_D"/env_vars.sh"
printf '#!/bin/sh\n\nunset PICARD\n' > $CONDA_DEACTIVATE_D"/env_vars.sh"

source deactivate

echo === The conda environment for crispr_scr is installed successfully. ===
