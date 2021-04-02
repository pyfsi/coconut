#!/bin/bash

# clean working directory
rm -rf ./CFD

# create new CFD folder
cp -r ../setup_files/openfoam3d CFD
cd CFD

# make blockmesh
ml OpenFOAM/4.1
source $FOAM_BASH
wmake ../../../coupling_components/solver_wrappers/openfoam/CoCoNuT_pimpleFoam/ > wmake_log
blockMesh > blockMesh_log

cd ..