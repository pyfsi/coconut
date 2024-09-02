#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
rm -rf ./CSM

# create new folder
mkdir CSM
cp -r setup_abaqus/nodes_elements.txt CSM
cp -r setup_abaqus/yarn.inp CSM