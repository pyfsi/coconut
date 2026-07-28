#!/bin/bash

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh_background.jou > setup_gambit_bg.log 2>&1
ml -GAMBIT

cp background_v2.msh ../
