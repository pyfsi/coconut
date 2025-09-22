#!/bin/bash

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh_bl.jou > setup_gambit.log 2>&1
ml -GAMBIT

cp 2D_liquid_BL.msh ../
