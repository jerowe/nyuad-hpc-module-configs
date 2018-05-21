#!/usr/bin/env bash

SYNC_DIR=$1
RUN=$2

aws s3 sync $1 `basename $1`
chmod 777 $RUN
bash $RUN
