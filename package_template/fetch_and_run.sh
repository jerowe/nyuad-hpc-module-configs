#!/usr/bin/env bash 

set -x -e

echo $USER
echo `pwd`
echo "hello from fetch and run..."

source ~/.bashrc

SYNC_DIR=$1
RUN=$2

mkdir -p hpc-runner/`basename $1`
aws s3 sync $1 hpc-runner/`basename $1`

chmod 777 $RUN
bash $RUN

#exec "$@"
