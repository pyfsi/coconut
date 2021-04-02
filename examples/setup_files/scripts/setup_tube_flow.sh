#!/bin/bash

solver=$1

# clean working directory
rm -rf ./CFD
mkdir ./CFD

# copy setup files
cp ../setup_files/tube_flow/solver_parameters_${solver}.json CFD/solver_parameters.json
