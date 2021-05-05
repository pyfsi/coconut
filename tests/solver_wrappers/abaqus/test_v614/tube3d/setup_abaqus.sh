#!/bin/bash

# README: run this script to remove old data and setup case

module load ABAQUS/6.14
module load intel/2018a

# clean working directory
rm -rf ./CSM

# create new folder
cp -r setup_abaqus CSM
cd CSM

source $PWD/../../../../../../coupling_components/solver_wrappers/abaqus/extra/make_host_file.sh

cd ..