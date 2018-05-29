#!/usr/bin/env bash 

source ~/.bashrc
source activate $EB_ENV
module load gencore_biosails

##REMOVE THIS BEFORE FINAL RELEASE
/home/ebuser/bin/update_hpc_runner.sh

set -x -e
echo "ARGS ARE"
echo "$@"

## This works on AWS
## But not on my local machine BAH
exec "${@}"

## Either of these works on my local machine but not on AWS
#bash -c "$@"
#FILE=$1
#REST_OF_ARGS="${@:2}"
#bash -c $FILE ${REST_OF_ARGS}
