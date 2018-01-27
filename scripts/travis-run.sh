#!/bin/bash

set -xe # Exit with nonzero exit code if anything fails

export BASE_DIR=`pwd`
# export PATH="$HOME/anaconda3/bin:$PATH"

git fetch origin master
#export RECIPES=$(git diff FETCH_HEAD --name-only | grep yml | grep recipes)
#echo "Processing recipes..."
#echo $RECIPES
#echo ""

##TRAVIS_BRANCH
##CIRCLE_BRANCH
# if [[ $TRAVIS_BRANCH = "master" && "$TRAVIS_PULL_REQUEST" = false ]]
# ? CIRCLE_COMPARE_URL
if [[ $CIRCLE_BRANCH = "master" ]]
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
    gencore_app build_envs
    # git diff FETCH_HEAD --name-only | grep yml | grep recipes | xargs -I {} gencore_app build_envs -e {}
fi

exit 0
