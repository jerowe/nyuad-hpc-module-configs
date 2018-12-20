#!/usr/bin/env bash

source /usr/share/lmod/lmod/init/bash
module load gencore_biosails

git clone -b feature/aws-batch https://github.com/biosails/HPC-Runner-Command
cd HPC-Runner-Command

sed -i.bak 's|#!/usr/bin/env perl|#!perl|g' script/hpcrunner.pl
rm -rf script/hpcrunner.pl.bak

HOME=/tmp cpanm --installdeps .
cpanm .

cd ..
rm -rf HPC-Runner-Command
