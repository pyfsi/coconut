# TestExamples

This documentation describes the different test examples.
Currently all these examples calculate the flow in a flexible tube.


## Folder and file structure

This section describes the different folders and files which are provided.

- run_simulation.py: Main file which has to be run with a parameterfile as argument
- project_parameters_X.json: Parameter file in json format
- setup_X: Setup folder containing all files for setting up solver X
- setup_X.sh: Bash file which has to be run to set up solver X
- X.md: A description of the specific example

When the setup files are run, working directories are created which have to match the ones specified in the parameter file.
These folder are expandable and are deleted when the setup files are (re)run.

## Project Parameters

The parameter file `project_parameters_X.json` contains two dictionaries. 
The first is `settings` in which the number of time steps for which the calculation is run is specified.
The second is `coupled_solver` in which all parameters related to solver algorithm are given: coupling algorithm, predictor, convergence criteria and solvers.

## Running a case

In order to run a test example, first, the setup files have to be run to create the working directories.
Then, the calculation is started by running run_simulation.py with the parameter file as argument.

## Debug files

The folder `test_examples` also contains a folder `debug_files` containing scripts for debug purposes.
These files might need some adjustements to work. 
In order to use them, the debug boolean `self.debug` has to be `True` in the code of the corresponding solver wrappers.