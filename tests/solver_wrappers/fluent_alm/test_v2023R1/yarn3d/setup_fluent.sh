#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./CFD ] 
then
    rm -rf ./CFD 
fi

# create new CFD folder
cp -r ./setup_fluent/ CFD
cd CFD
gunzip mesh_yarn3d.msh.gz

# make fluent case
fluent 3ddp -g -i case.jou > setup_fluent.log 2>&1
