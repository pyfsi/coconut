#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory
cp ../setup_files/run_simulation.py ./

# clean working directories
rm -rf ./CFD
rm -rf ./CSM

# create new CFD folder
cp -r ../setup_files/tube/tube_flow CFD

# create new CSM folder
cp -r ../setup_files/tube/abaqus2d CSM
cd CSM; ./setup_abaqus2d.sh; cd ..

# TODO: remove once each solver runs in own terminal
module load ABAQUS/6.14
module load intel/2018a
