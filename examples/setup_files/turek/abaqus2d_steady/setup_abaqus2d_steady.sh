#!/bin/bash

# create model
# TODO remove module load command
module load ABAQUS/6.14
abaqus cae noGUI=make_inp.py
rm *.rpy*

source $PWD/../../../coupling_components/solver_wrappers/abaqus/extra/make_host_file.sh

cd ..
