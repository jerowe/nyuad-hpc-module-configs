#!/usr/bin/env bash

set -x -e

source activate $EB_ENV
module load gencore_biosails

git clone -b feature/aws-batch https://github.com/biosails/HPC-Runner-Command
cd HPC-Runner-Command
cpanm .
