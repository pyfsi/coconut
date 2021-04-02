#!/bin/bash

# clean working directory
rm -rf ./CFD

# create new CFD folder
cp -r ../setup_files/fluent2d CFD
cd CFD

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou > setup_gambit.log 2>&1
ml -GAMBIT

# make fluent case
ml ANSYS_CFD/2019R1
fluent 2ddp -g -i case.jou > setup_fluent.log 2>&1
ml -ANSYS_CFD

cd ..

ml ANSYS_CFD/2019R1
