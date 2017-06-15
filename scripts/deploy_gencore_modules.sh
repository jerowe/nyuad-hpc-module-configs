#!/usr/bin/env bash

### NYUAD CGSB EasyBuild Module Deployment
##################################################################
## Base Image Setup
##################################################################

#sudo yum install -y git wget tar make bzip2 patch which
#sudo apt-get install -y git wget tar make bzip2 patch which

##################################################################
## Anaconda Setup
## We are using miniconda3/python3 from anaconda
## Be sure to change this if you are running as a your own user
##################################################################
export ANACONDA3_BASE=/home/jillian/.local/anaconda3

curl -s -L https://repo.continuum.io/miniconda/Miniconda3-4.2.12-Linux-x86_64.sh  > miniconda.sh && \
    openssl md5 miniconda.sh | grep d0c7c71cc5659e54ab51f2005a8d96f3

chmod 777 miniconda.sh ; bash miniconda.sh -b -p $ANACONDA3_BASE 
rm miniconda.sh && \
    export PATH=${ANACONDA3_BASE}/bin:$PATH && \
    conda config --set show_channel_urls True && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda update --all --yes && \
    conda clean -tipy

##################################################################
## EasyBuild User Setup
## This is only necessary if you are creating a user for EB
## If running as your own user skip this section
## But you cannot run easybuild as root
##################################################################
#useradd -ms /bin/bash ebuser
#
#cd /home/ebuser
#
#export HOME=/home/ebuser
#
##################################################################
## EasyBuild User Setup
## This is only necessary if you
##################################################################

export PATH=${ANACONDA3_BASE}/bin:$PATH
export EASYBUILD_PREFIX=$HOME/.eb
export EB_ENV=$HOME/.eb/eb--3.1.10_lmod--7.3.16
export MODULEPATH=$EASYBUILD_PREFIX/modules/all

##################################################################
## Add Easybuild Config
## For now the config doesn't have anything except some commented
## out sections for the modules-tool and modules-syntax
##################################################################
mkdir -p $HOME/.config/easybuild

##################################################################
## Conda config for ebuser
## Setup the initial conda configuration for the ebuser
##################################################################
conda config --set show_channel_urls True && \
    conda config --add channels nyuad-cgsb && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels r && \
    conda config --add channels bioconda && \
    conda config --set always_yes True && \
    conda config --set allow_softlinks False


##################################################################
## Install Easybuild and Lmod
## Create a conda env with easybuild and lmod - and nothing else
##################################################################
conda create -p $EB_ENV easybuild lmod

##################################################################
## Install an Easybuild Config
## Anaconda2 is a Core EB config
## --software-name is a way of telling EB to look for software
## from the repos instead of supplying a .eb file
##################################################################
source activate $EB_ENV && \
	eb --software-name Anaconda2

##################################################################
## Add NYUAD Easybuild Configs
## These are custom configs outside of the Easybuild Main
##################################################################
mkdir -p $HOME/.eb/custom_repos
cd $HOME/.eb/custom_repos
git clone https://github.com/nyuad-cgsb/nyuad-hpc-module-configs
cd $HOME

###################################################################
### Add NYUAD Easybuild Configs to the Robot Path
### Robot will tell EB to automatically pull in deps
### Robot-path will tell EB where to look for configs
###################################################################
export ROBOT=$HOME/.eb/custom_repos/nyuad-hpc-module-configs/_easybuild

### Remove the --dry-run in order to actually install the module
source activate $EB_ENV
eb --dry-run  --robot --robot-paths=$ROBOT   gencore_rnaseq-1.0.eb

### Using --extended-dry-run will give you more information
eb --extended-dry-run  --robot --robot-paths=$ROBOT  gencore_rnaseq-1.0.eb

module avail

##################################################################
## Install a gnecore module
## Finally, running this command will install the gencore_rnaseq-1.0 module
## Additionally, it will install the gencore_anaconda2 module if it is not already available
##################################################################

# eb  --robot --robot-paths=$ROBOT  gencore_rnaseq-1.0.eb
# module avail
