#!/bin/bash

# add coconut installation path to parent directory of coconut
COCONUT_INSTALLATION_PATH=$(dirname $(pwd))
export PYTHONPATH=$COCONUT_INSTALLATION_PATH:$PYTHONPATH

#run tests
echo "Running tests..."
cd $COCONUT_INSTALLATION_PATH/coconut/tests
python3 test_coconut.py

# run an FSI test example
echo "Running test example..."
cd $COCONUT_INSTALLATION_PATH/coconut/test_examples/tube_tube_flow_tube_structure
chmod 700 setup_tube_flow.sh setup_tube_structure.sh
./setup_tube_flow.sh
./setup_tube_structure.sh
python3 run_simulation.py project_parameters_mapped.json > run.log

# installation successful
echo "CoCoNuT installation successful."