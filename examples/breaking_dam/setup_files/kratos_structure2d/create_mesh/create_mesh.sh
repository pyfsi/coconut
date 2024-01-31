#!/bin/bash

source ~/.bashrc

addtool gmsh
gmsh beam_mesh.geo -2 -format msh2 -o beam_mesh.msh

ml Anaconda3-python
python gmsh_to_mdpa.py beam_mesh.msh mesh_parameters.json
ml - Anaconda3-python

cp beam_mesh.mdpa ../mesh_breaking_dam.mdpa
