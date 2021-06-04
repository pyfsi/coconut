#!/bin/bash

# create model
abaqus cae noGUI=make_inp.py > setup_abaqus.log 2>&1
rm *.rpy*
