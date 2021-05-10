#!/bin/bash

# make gambit mesh
#TODO remove module load command
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh.jou > setup_gambit.log 2>&1
ml -GAMBIT
