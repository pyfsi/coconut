#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
rm -rf ./CSM

# create new folder
cp -r ../../test_v614/tube2d/setup_abaqus CSM
