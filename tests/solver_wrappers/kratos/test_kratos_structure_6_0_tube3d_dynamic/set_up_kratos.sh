#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./kratos_structure ]
then
    rm -rf ./kratos_structure
fi

# create new CSM folder
cp -r kratos_setp_up_files kratos_structure
cd kratos_structure

ml Kratos
cd ..