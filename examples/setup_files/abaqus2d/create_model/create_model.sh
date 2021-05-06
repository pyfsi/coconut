#!/bin/sh 

abaqus cae noGUI=makeInp.py

mv CSM_Time0.inp ../

rm *.rpy*
