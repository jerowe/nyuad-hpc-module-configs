#!/bin/bash

# set -euo pipefail
#
# if [[ $TRAVIS_OS_NAME = "linux" ]]
# then
#
#     docker run -e TRAVIS_PULL_REQUEST -e TRAVIS_BRANCH \
#         -e GITHUB_TOKEN -e GITHUB_USER -e GITHUB_REPO \
#         -e ANACONDA_TOKEN -e ANACONDA_PASSWORD -e ANACONDA_USER \
#         -i -t -v `pwd`:/nyuad-conda-configs quay.io/nyuad_cgsb/anaconda-centos /nyuad-conda-configs/tests/run-tests.sh
#
# else
#     exit 0
# fi
set -xe # Exit with nonzero exit code if anything fails

export BASE_DIR=`pwd`

export PATH="/home/travis/anaconda3/bin:$PATH"
/home/travis/anaconda3/bin/pip install git+https://github.com/nyuad-cgsb/gencore_app.git@master

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
    # cd /nyuad-conda-configs
    cd $BASE_DIR
    gencore_app build_eb

    echo "Gencore App building docs"
    # cd /nyuad-conda-configs
    cd $BASE_DIR
    gencore_app build_docs

    echo "Commit docs to github"
    # cd /nyuad-conda-configs
    cd $BASE_DIR
    scripts/build_easybuild.sh
    scripts/build_docs.sh

    echo "Building man pages!"
    # cd /nyuad-conda-configs
    cd $BASE_DIR
    gencore_app build_man

    echo "Uploading packages to anaconda!"
    gencore_app upload_envs

    scripts/commit_recipes.sh
else
    #Just test packages
    gencore_app --help
    gencore_app build_envs
    # git diff FETCH_HEAD --name-only | grep yml | grep recipes | xargs -I {} gencore_app build_envs -e {}
fi
