#!/bin/bash

# remesh
ml -ANSYS_CFD
ml GAMBIT/2.4.6
gambit -inp remesh_tri.jou > setup_gambit.log 2>&1
ml -GAMBIT
