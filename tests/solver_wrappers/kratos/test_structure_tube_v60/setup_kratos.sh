#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./CSM ]
then
    rm -rf ./CSM
fi

# create new CSM folder
cp -r setup_kratos CSM

