#!/bin/bash

# clean working directory
rm -rf ./CSM
mkdir ./CSM

solver=$1

# copy setup files
cp ../setup_files/tube_structure/solver_parameters_${solver}.json CSM/solver_parameters.json