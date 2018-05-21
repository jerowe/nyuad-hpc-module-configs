#!/usr/bin/env bash

set -x -e

echo `pwd`
echo "hello from fetch and run..."

SYNC_DIR=$1
RUN=$2

mkdir -p `basename $1`
aws s3 sync $1 hpc-runner/`basename $1`
chmod 777 $RUN
bash $RUN

ls -lah `basename $1`
