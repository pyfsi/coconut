#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory
cp ../setup_files/run_simulation.py ./

# clean working directories
rm -rf ./CFD
rm -rf ./CSM

# create new CFD folder
cp -r ../setup_files/openfoam3d CFD
cd CFD; ./setup_openfoam3d.sh; cd ..

# create new CSM folder
cp -r ../setup_files/kratos_structure3d CSM
