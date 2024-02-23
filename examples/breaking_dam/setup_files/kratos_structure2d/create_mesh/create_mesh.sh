#!/bin/bash

ml gmsh
gmsh beam_mesh.geo -2 -format msh2 -o beam_mesh.msh
ml purge

ml Anaconda3-python
 # gmsh_to_mdpa.py should be installed and found in PYTHONPATH, see https://github.com/pyfsi/gmsh_to_mdpa
python gmsh_to_mdpa.py beam_mesh.msh mesh_parameters.json
ml - Anaconda3-python

cp beam_mesh.mdpa ../mesh_breaking_dam.mdpa
