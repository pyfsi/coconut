#!/bin/bash

# make gambit mesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp mesh_background.jou > setup_gambit_bg.log 2>&1
gambit -inp mesh_component.jou > setup_gambit_co.log 2>&1
ml -GAMBIT

cp background.msh ../
cp component.msh ../
