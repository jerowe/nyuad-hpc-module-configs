#!/usr/bin/env bash 

source ~/.bashrc
source /usr/share/lmod/lmod/init/bash
source activate $EB_ENV
module load gencore_biosails

##REMOVE THIS BEFORE FINAL RELEASE
/home/ebuser/bin/update_hpc_runner.sh

set -x -e

## This works on AWS
## But not on my local machine BAH
if [ -z ${AWS_BATCH_JOB_ID+x} ]; then
	echo "This is an AWS Batch job..."
	exec "${@}"
else 
	echo "var is set to '$var'"; 
	bash -c "$@"
fi

## Either of these works on my local machine but not on AWS
#FILE=$1
#REST_OF_ARGS="${@:2}"
#bash -c $FILE ${REST_OF_ARGS}
