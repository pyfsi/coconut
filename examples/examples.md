# Examples

This documentation describes the structure of the different example cases and how to run them.


## Folder and file structure

The *`examples`* folder contains different examples grouped per case, for example, the cases calculating the flow through a flexible tube is gathered in the folder *`tube`*.

Each case folder, such as *`tube`*, contains the following directories and files:

- *`setup_files`* containing the necessary files to set up the case for each solver,
- *`<flow solver>_<structural solver>`* different folders solving the case with different solver combinations, hereafter simply called examples.

The following files can be found in each example:

- *`parameters.json`* the parameter file in JSON format specifying the settings for the CoCoNuT simulation, and which can be interpreted as a dictionary (there are two main keys,
the first is `settings` in which the number of time steps for which the calculation is run is specified and the second is `coupled_solver` in which all parameters related to solver algorithm are given: coupling algorithm, predictor, convergence criteria and solver wrappers),
- *`setup_case.py`* the Python script used to create the working directories (copying them from *`setup_files`*, mentioned earlier), to set up the case for the different solvers and to copy *`run_simulation.py`* (note that these folders are expandable and are overwritten when the setup file is (re)run), 
- *`tube_<flow solver>_<structural solver>.md`* description of the specific example,
- optionally other files which can be used for visualization of the interface data or for comparison with benchmark data.

Besides the case folders, such as *`tube`*, the *`examples`* folder contains other directories:

- *`test_single_solver`* particular example of the coupled solver `test_single_solver` that can be used to test the setup of case by only running one solver and replacing the other solver with a dummy solver that outputs prescribed data,
- *`post_processing`* containing the code to inspect and visualize your results, see also the [post-processing documentation](post_processing/post_processing.md), given that your results are stored in a pickle file, see the [documentation on save results](../coupling_components/coupled_solvers/coupled_solvers.md#save-results).
- *`run_simulation.py`* the script used to start the actual simulation,
- *`images`* folder containing images used in the description of the examples,
- *`evaluate_examples`* Python script to compare the results of the example simulation with benchmark data,
- *`benchmark_files`* folder containing the benchmark data in the form of pickle files and a Python script to create this data.


## Running a case

In order to run an example, first, navigate to the example and run the setup script to create the working directories and set up the cases for each solver.
This script also copies *`run_simulation.py`* to the example.
The calculation is started by running *`run_simulation.py`*.

In practice, you will run the following commands:
```bash
python3 setup_case.py
python3 run_simulation.py
```

**Important**
CoCoNuT launches each solver in its own environment, avoiding potential conflicts.
This means that CoCoNuT needs to know the commands to make the software available on your machine.
These commands are specified in the *`coconut/solver_modules.py`* file, in the dictionary `solver_load_cmd_dict`.
It's important to add your own machine to this dictionary and to specify the `machine_name` correspondingly.
The solvers for the UGent-cluster have already been added to the file.
If it's not necessary to, for example, load modules, on your machines, an empty string should be used.


## Setting up your own case

To set up your own case, following items should be present:

- a JSON file containing the parameters of your simulation; in the examples, this file is called *`parameters.json`*, but you can deviate from this convention as long as the variable `parameter_file_name` in the *`run_simulation.py`* script is modified,
- a *`run_simulation.py`* script which starts the CoCoNuT simulation when run,
- working directories for the used solvers, contain the case files required for the calculations; in the examples they are named *`CFD`* and *`CSM`*, respectively for the flow and structural solvers, but they can have any name as long as they are referenced correctly in the parameter file under the keyword `working_directory`.

Finally, some tips are listed to successfully set up you and run your own case:

- Although it is not required, it can be a good idea to set up your case using a set-up-script and folder with setup-files as in the examples, allowing you to easily and reproducibly start from a clean case.
- The easiest way to get started is to start from an existing example closest to your own case and modify the files as required.
- Before coupling your solvers, it is wise to test them first separately using `test single solver` and a `dummy solver`, which can prescribe, for example, displacement of the interface. This way it can be verified if the solvers behave as expected before coupling them together.
- The scripts for animating the data on the interface provided in *`post_processing`* can be a handy tool to understand possible issues.
