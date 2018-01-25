#!/usr/bin/env bash

set -ex

echo "BUILDING EASYBUILD"

if [[ -z "${GITHUB_TOKEN}" ]] ; then
    echo "GitHub API key needs to be set to update docs."
    exit 0
fi

#cd /nyuad-conda-configs

#At least we can test if this works
##TODO Make these separate github repos
mkdir -p _easybuild

git status
git add -A

git config  user.email "jillian.e.rowe@gmail.com"
git config  user.name "jerowe"

ORIGIN="https://${GITHUB_USER}:${GITHUB_TOKEN}@github.com/${GITHUB_USER}/${GITHUB_REPO}.git"

git remote rm origin
git remote add origin "$ORIGIN"

#git add -A
# git pull origin "$TRAVIS_BRANCH"
git pull origin "$CIRCLE_BRANCH"

git add _easybuild
#If it doesn't exit as 0 there is nothing to commit
# git commit  -m "Updated docs to commit ${TRAVIS_COMMIT}." || exit 0
git commit  -m "Updated docs to commit" || exit 0
# git push  origin "$TRAVIS_BRANCH"
git push  origin "$CIRCLE_BRANCH"
