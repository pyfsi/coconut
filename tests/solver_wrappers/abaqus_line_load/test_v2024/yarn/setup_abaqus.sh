#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
rm -rf ./CSM

# create new folder
mkdir CSM
cp -r ../../test_v2024/yarn/setup_abaqus/nodes_elements.txt CSM
cp -r ../../test_v2024/yarn/setup_abaqus/yarn.inp CSM