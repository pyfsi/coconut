#!/bin/sh

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou > setup_gambit.log 2>&1
ml -GAMBIT

python3 convert_neu_to_inp.py mesh_breaking_dam.neu mesh_breaking_dam.inp CPE8R

cp mesh_breaking_dam.inp ../
