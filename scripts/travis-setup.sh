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

sudo chown -R $USER /opt/anaconda3

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh --quiet -O miniconda.sh
bash miniconda.sh -b -p /home/travis/anaconda3

export PATH="/home/travis/anaconda3/bin:$PATH"

conda config --set always_yes yes --set changeps1 no
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nyuad-cgsb

# conda index /home/travis/anaconda3/conda-bld/linux-64
# conda config --add channels file://home/travis/anaconda3/conda-bld
conda update --all -y

conda install python=3.5
conda install conda conda-build
conda install -y r-base r-essentials openjdk perl bioconductor-biobase nodejs
npm install -g marked-man

conda install -y pip
/home/travis/anaconda3/bin/pip install git+https://github.com/nyuad-cgsb/gencore_app.git@master
