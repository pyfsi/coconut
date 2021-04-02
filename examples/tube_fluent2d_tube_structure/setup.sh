#!/bin/bash

# README: run this script to remove old data and setup case

# copy run_simulation.py script to main directory

cp ../setup_files/run_simulation.py ./

# setup CSM folder
source ../setup_files/scripts/setup_tube_structure.sh fluent

# setup CFD folder
source ../setup_files/scripts/setup_fluent2d.sh