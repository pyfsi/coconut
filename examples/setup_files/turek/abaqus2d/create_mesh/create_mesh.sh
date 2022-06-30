#!/bin/sh

#TODO remove module load command
# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou > setup_gambit.log 2>&1
ml -GAMBIT

python convert_neu_to_inp.py mesh_turek.neu mesh_turek.inp CPE8R

cp mesh_turek.inp ../
