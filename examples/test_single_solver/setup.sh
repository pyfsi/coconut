#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory
cp ../setup_files/run_simulation.py ./

# setup CSM
source ../setup_files/scripts/setup_abaqus3d.sh

# setup CFD
source ../setup_files/scripts/setup_fluent3d.sh