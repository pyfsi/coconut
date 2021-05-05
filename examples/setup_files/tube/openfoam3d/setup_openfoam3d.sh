#!/bin/bash

#TODO remove module load command
ml OpenFOAM/4.1
source $FOAM_BASH
wmake ../../../coupling_components/solver_wrappers/openfoam/CoCoNuT_pimpleFoam/ > wmake_log
blockMesh > blockMesh_log
cd ..
