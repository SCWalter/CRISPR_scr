#!/bin/bash
# Stop on error
set -e

CUR_DIR=$(pwd)

conda create -n fast_iclip --file requirements.txt -y -c defaults -c bioconda -c r

############ install additional packages

source activate fast_iclip

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
mkdir -p $CONDA_EXTRA

#### install htseq==0.6.1
pip install htseq==0.6.1

#### install setuptools_git
cd $CONDA_EXTRA
wget https://github.com/msabramo/setuptools-git/archive/1.1.tar.gz
tar -xvzf 1.1.tar.gz
rm -f 1.1.tar.gz
cd setuptools-git-1.1
python setup.py install --record files.txt

#### install CLIPper
cd $CONDA_EXTRA
wget https://github.com/YeoLab/clipper/archive/1.1.tar.gz
tar -xvzf 1.1.tar.gz
rm -f 1.1.tar.gz
cd clipper-1.1
python setup.py install --record files.txt

#### install iCLIPro

cd $CONDA_EXTRA
wget http://www.biolab.si/iCLIPro/dist/iCLIPro-0.1.1.tar.gz -N
tar -xvzf iCLIPro-0.1.1.tar.gz
rm -f iCLIPro-0.1.1.tar.gz
cd iCLIPro-0.1.1
python setup.py install --record files.txt

cd $CUR_DIR

# 2. Downloading genome files

echo "Downloading and unzipping hg19 and mm9 genome files"
wget https://s3.amazonaws.com/changlabguest/FAST-iCLIP/docs.tar.gz
mkdir docs
tar xvzf docs.tar.gz -C docs && rm -f docs.tar.gz

# 3. Downloading example files
wget https://s3.amazonaws.com/changlabguest/FAST-iCLIP/rawdata.tar.gz
mkdir rawdata
tar xvzf rawdata.tar.gz -C rawdata && rm -f rawdata.tar.gz
mkdir results

echo === Installing dependencies successfully done. ===

python setup.py install --record files.txt

source deactivate

echo === Installing FAST-iCLIP successfully done. ===
