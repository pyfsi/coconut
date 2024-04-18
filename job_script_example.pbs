#!/bin/sh
#
#PBS -N example_lid_driven_cavity
#PBS -l walltime=00:59:00
#PBS -l nodes=1:ppn=all
#

ml Anaconda3
export PYTHONPATH=$HOME/Tools/coconut/:$PYTHONPATH
export PYTHONPATH=$HOME/Tools/Kratos/lib/python3.10/site-packages/:$PYTHONPATH

export OMP_PROC_BIND=FALSE # required on UGent-HPC for processes spawned from Python

cd $HOME/Tools/coconut/coconut/examples/lid_driven_cavity/fluent2d_kratos_structure2d/
python3 setup_case.py
python3 run_simulation.py

