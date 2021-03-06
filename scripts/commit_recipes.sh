#!/usr/bin/env bash

set -ex

if [[ -z "${GITHUB_TOKEN}" ]] ; then
    echo "GitHub API key needs to be set to update docs."
    exit 0
fi

#cd /nyuad-conda-configs

#At least we can test if this works
echo "Recommiting recipes"

git status

git config  user.email "jillian.e.rowe@gmail.com"
git config  user.name "jerowe"

git add _docs

ORIGIN="https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/${GITHUB_USER}/${GITHUB_REPO}.git"
git remote rm origin
git remote add origin "$ORIGIN"

git add recipes
#IF it doesn't exit as 0 its because there is nothing to commit
git commit  -m "Updated recipes " || exit 0
# git push -f origin "$TRAVIS_BRANCH"
git push -f origin "$CIRCLE_BRANCH"
