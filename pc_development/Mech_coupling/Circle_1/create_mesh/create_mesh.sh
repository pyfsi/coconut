#!/bin/bash

# prepare gambit script (to run on cfdclu13)
ml Anaconda3-python
python3 initiate.py

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh_filled.jou > setup_gambit.log 2>&1
ml -GAMBIT

cp liquid.msh ../
