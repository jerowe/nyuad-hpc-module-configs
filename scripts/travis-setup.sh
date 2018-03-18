#!/bin/bash
set -e

# if [[ $TRAVIS_OS_NAME = "linux" ]]
# then
#     docker pull quay.io/nyuad_cgsb/anaconda-centos
# else
#
#     exit 0
# fi

# pip install pyyaml

# sudo chown -R $USER /opt/anaconda3

# wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --quiet -O miniconda.sh
# bash miniconda.sh -b -p $HOME/anaconda3
#
# export PATH="$HOME/anaconda3/bin:$PATH"

apt-get update -y
apt-get install -y  build-essential

conda config --set always_yes yes --set changeps1 no
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nyuad-cgsb

mkdir -p $HOME/anaconda3/conda-bld

conda update --all -y

##These are all added to the travis cache

conda install python=3.5
conda install -y conda conda-build anaconda-client pip setuptools
conda install -y r-base nodejs 
conda install -y gnuplot samtools bamtools bcftools freebayes gatk
conda install -y openjdk perl bioconductor-biobase blast bedtools
npm install -g marked-man

pip uninstall gencore_app || echo "Gencore app is not installed"
pip install git+https://github.com/nyuad-cgsb/gencore_app.git@master
