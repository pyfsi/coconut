#!/bin/bash

export INTEL_LICENSE_FILE=28518@157.193.126.6
source /apps/SL6.3/Intel/compiler/2015.3.187/bin/compilervars.sh intel64
module load ABAQUS/6.14

# clean working directory
rm -rf ./CSM

# create new folder
cp -r ../setup_files/abaqus2d_steady CSM
cd CSM

source makeHostFile.sh

cd ..