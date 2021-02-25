#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./CSM ]
then
    echo "Removing kratos_structure folder"
    rm -rf ./CSM
fi

# create new CSM folder
echo "Copying kratos_structure and set up files"
cp -r setup_kratos/ CSM
