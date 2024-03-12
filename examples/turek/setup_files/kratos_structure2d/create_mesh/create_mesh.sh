#!/bin/bash

ml gmsh
gmsh flag_mesh.geo -1 -2 -format msh2 -o flag_mesh.msh
ml purge

ml Anaconda3-python
 # gmsh_to_mdpa.py should be installed and found in PYTHONPATH, see https://github.com/pyfsi/gmsh_to_mdpa
gmsh_to_mdpa flag_mesh.msh mesh_parameters.json
ml -Anaconda3-python

cp flag_mesh.mdpa ../flag.mdpa
