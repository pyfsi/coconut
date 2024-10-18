#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./CFD ] 
then
    rm -rf ./CFD 
fi

# create new CFD folder
cp -r ../../test_v2024R1/tube2d/setup_fluent/ CFD
cd CFD

# make fluent case
fluent 2ddp -g -i case.jou > setup_fluent.log 2>&1
