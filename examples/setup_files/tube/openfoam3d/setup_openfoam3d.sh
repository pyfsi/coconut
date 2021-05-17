#!/bin/bash

#compile CoCoNuT_pimpleFoam
wmake ../../../coupling_components/solver_wrappers/openfoam/coconut_pimpleFoam/ > wmake_log

#create mesh
blockMesh > blockMesh_log
