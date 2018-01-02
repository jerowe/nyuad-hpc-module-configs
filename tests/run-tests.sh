#!/bin/bash

set -xe # Exit with nonzero exit code if anything fails

export PATH=/anaconda/bin:$PATH

conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda


conda update -y

conda install -y pip
conda install conda conda-build


pip install git+https://github.com/nyuad-cgsb/gencore_app.git@master


cd /nyuad-conda-configs

#git fetch origin master
#export RECIPES=$(git diff FETCH_HEAD --name-only | grep yml | grep recipes)
#
#echo "Processing recipes..."
#echo $RECIPES
#echo ""

if [[ $TRAVIS_BRANCH = "master" && "$TRAVIS_PULL_REQUEST" = false ]]
then
    #Upload packages
    anaconda login --user $ANACONDA_USER --password $ANACONDA_PASSWORD
    conda config --set anaconda_upload yes

    echo "Gencore App build ebs"
    cd /nyuad-conda-configs
    gencore_app build_eb

    echo "Gencore App building docs"
    cd /nyuad-conda-configs
    gencore_app build_docs

    echo "Commit docs to github"
    cd /nyuad-conda-configs
    scripts/build_easybuild.sh
    scripts/build_docs.sh

    echo "Building man pages!"
    cd /nyuad-conda-configs
    gencore_app build_man

    echo "Uploading packages to anaconda!"
    gencore_app upload_envs

    scripts/commit_recipes.sh
else
    #Just test packages
    gencore_app build_envs
fi
