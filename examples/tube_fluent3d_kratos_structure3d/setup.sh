#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory
cp ../setup_files/run_simulation.py ./

# clean working directories
rm -rf ./CFD
rm -rf ./CSM

# create new CFD folder
cp -r ../setup_files/tube/fluent3d CFD
cd CFD; ./setup_fluent3d.sh; cd ..

# create new CSM folder
cp -r ../setup_files/tube/kratos_structure3d CSM

# TODO: remove once each solver runs in own terminal
module load ANSYS_CFD/2019R1
