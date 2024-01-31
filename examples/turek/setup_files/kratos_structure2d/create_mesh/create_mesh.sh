#!/bin/bash

source ~/.bashrc

addtool gmsh
gmsh flag_mesh.geo -1 -2 -format msh2 -o flag_mesh.msh

ml Anaconda3-python
addpath gmsh_to_mdpa
gmsh_to_mdpa flag_mesh.msh mesh_parameters.json
ml -Anaconda3-python

cp flag_mesh.mdpa ../flag.mdpa
