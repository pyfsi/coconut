#!/bin/bash

# README: run this script to remove old data and setup case

# clean working directory
if [ -d ./kratos_structure ]
then
    echo "Removing kratos_structure folder"
    rm -rf ./kratos_structure
fi
rm *.msh
rm *.res
rm *.lst
rm *.time
# create new CSM folder
echo "Copying kratos_structure and set up files"
cp -r kratos_setp_up_files kratos_structure
cd kratos_structure

cd ..
