#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory
cp ../setup_files/run_simulation.py ./

# clean working directories
rm -rf ./CFD
rm -rf ./CSM

# create new CFD folder
cp -r ../setup_files/fluent3d CFD
cd CFD; ./setup_fluent3d.sh; cd ..

# create new CSM folder
cp -r ../setup_files/abaqus2d CSM
cd CSM; source makeHostFile.sh; cd ..

# TODO: remove once each solver runs in own terminal
module load ANSYS_CFD/2019R1
module load ABAQUS/6.14
module load intel/2018a
