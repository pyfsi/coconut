#!/bin/bash

# README: run this script to setup case

# clean working directory
rm -rf ./CFD
mkdir ./CFD

# copy setup files
cp setup_tube_flow/solver_parameters.json CFD/