#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./CFD ] 
then
    rm -rf ./CFD 
fi

# create new CFD folder
cp -r ./../../test_v2019R3/tube3d/setup_fluent/ CFD
cp ./setup_fluent/case.jou CFD  # overwrite case.jou file
cd CFD

# make fluent case
fluent 3ddp -g -i case.jou > setup_fluent.log 2>&1
