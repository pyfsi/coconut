# Examples

This documentation describes the different example cases.
Currently all these examples calculate the flow in a flexible tube.


## Folder and file structure

This section describes the different folders and files that are provided.

The *`examples`* folder contains following subfolders:

- *`debug_files`*: Python scripts that are useful for debugging of the cases. These files might need some adjustments to work. In order to use them, the debug boolean `self.debug` has to be set to `True` in the settings of the corresponding solver wrappers.
- *`post_processing`*: contains scripts to visualize your results. In particular, the *`animate_example.py`* script provides the means to animate the variables stored in the interface (pressure, traction, displacement), given that your results are stored in a pickle file, see the [documentation on save results](../coupling_components/coupled_solvers/coupled_solvers.md#save-results).
- *`setup_files`*: contains the necessary files to set up the case for each solver and the *`run_simulation.py`* script necessary to start the actual simulation.
- *`test_single_solver`*: example to check the correct set up of your solver specific cases.
- *`tube_<flow solver>_<structural solver>`*: main directory of the example cases.

The main directory of each example case contains following files: 

- *`parameters.json`*: parameter file in JSON format.
- *`setup_case.py`*: Python script, which needs to be run to set up the case.
- *`tube_<flow solver>_<structural solver>.md`*: description of the specific example.

When the setup file is run, working directories are created that have to match the ones specified in the parameter file.
These folders are expandable and are deleted when the setup files are (re)run. 
The script copies first the *`run_simulation.py`* file into the main directory and then the necessary case files from *`setup_files`* into the correct working directory.
After the completion of this script, the example case is ready to run.


## Project Parameters

The parameter file *`parameters.json`* can be interpreted as a dictionary. 
There are two main keys.
The first is `settings` in which the number of time steps for which the calculation is run is specified.
The second is `coupled_solver` in which all parameters related to solver algorithm are given: coupling algorithm, predictor, convergence criteria and solvers.


## Running a case

In order to run an example case, first, the setup script has to be run to create the working directories.
Then, the calculation is started by running *`run_simulation.py`*. Make sure to add the module-load commands for the required solvers in `solverload_cmd_dict` in the *`coconut/solver_modules.py`* file, to load the necessary modules for your solvers during the simulation. By default, the solvers for the UGent-cluster are added to the file.

In practice, you would run the following commands:

```bash
python3 setup_case.py
python3 run_simulation.py
```

## Setting up your own case

To set up your own case, following items should be present in the main directory:

- A JSON file containing the parameters of your simulation. In the example cases, it is called *`parameters.json`*. 
In case you wish another name, make sure to adapt the variable `parameter_file_name` in the *`run_simulation.py`* script.
- The *`run_simulation.py`* script which initiates the CoCoNuT simulation when run.
- Working directories for the flow and structural solvers. In the examples they are named *`CFD`* and *`CSM`*, respectively. 
In practice, they can bear any name as long as you make sure they are referenced correctly in the parameter file under the keyword `working_directory`.
These working directories should contain the case files for the flow and structural calculations.