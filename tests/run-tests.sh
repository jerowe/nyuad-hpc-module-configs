#!/bin/bash

set -e # Exit with nonzero exit code if anything fails

export PATH=/anaconda/bin:$PATH

cd /nyuad-conda-configs


if [[ $TRAVIS_BRANCH = "master" && "$TRAVIS_PULL_REQUEST" = false ]]
then
    #Upload packages
    #TODO One script for testing and one for uploading

    anaconda login --user $ANACONDA_USER --password $ANACONDA_PASSWORD
    echo "Uploading packages to anaconda!"
    python3 tests/test_environments.py --master

    #echo "Uploading packages to anaconda!"
    #python3 tests/test_environments.py --master --environment recipes/variant_detection/1.0/environment-1.0.yml --force_rebuild

    #We will add this back soon
    #python3 scripts/build_docs.py --master
    #scripts/build_docs.sh
else
    #Just test packages
    python3 tests/test_environments.py

    #Protocol for forcing a package and verbose
    #python3 tests/test_environments.py --force_rebuild --verbose --environments recipes/qc/1.0/environment-1.0.yml
fi
