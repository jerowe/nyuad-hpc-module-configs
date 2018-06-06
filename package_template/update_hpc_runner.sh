#!/usr/bin/env bash

source /usr/share/lmod/lmod/init/bash
module load gencore_biosails

git clone -b feature/aws-batch https://github.com/biosails/HPC-Runner-Command
cd HPC-Runner-Command
cpanm .

cd ..
rm -rf HPC-Runner-Command
