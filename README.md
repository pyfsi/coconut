[![CoCoNuT banner](https://raw.githubusercontent.com/pyfsi/coconut/master/docs/images/coconut-banner.svg)](https://github.com/pyfsi/coconut)


# Coupling Code for Numerical Tools


CoCoNuT is a light-weight Python package for efficient partitioned multi-physics simulations, with a focus on fluid-structure interaction. 
Thanks to its fully modular approach, the package is versatile and easy to extend. It is available under the GPL-3.0 license. 


## Introduction

The *Coupling Code for Numerical Tools* &mdash; *CoCoNuT* in short &mdash; follows a partitioned approach to solve multi-physics problems: existing single-physics solvers are coupled through a Python interface. 
This has the advantage that dedicated, highly-optimized single-physics solvers can be used. 
To use the code with a new solver (open source or commercial), a so-called *solver wrapper* is written to take care of the communication with CoCoNuT. 
All other CoCoNuT components, such as mapping for non-conformal meshes, are solver-independent and are easily swapped thanks to CoCoNuT's modular architecture. 

CoCoNuT is under active development by the [Fluid Mechanics research team](https://www.ugent.be/ea/eemmecs/en/research/stfes/flow) at Ghent University. 
Our specialization is partitioned fluid-structure interaction. 
We develop high-performance quasi-Newton algorithms to accelerate convergence of the coupled problem, and apply these techniques to diverse fluid-structure interaction applications such as wind turbines, tube bundles and flexible aircraft. 

The full documentation of the CoCoNuT package can be found at the [documentation website](https://pyfsi.github.io/coconut/).


## Installation

These instructions describe the setup of CoCoNuT on Linux. The package has not been tested on Windows or macOS, so compatibility is not guaranteed, although we do not expect major issues. 


### Requirements

-   `python>=3.11.5` 
-   `numpy>=1.24.3`
-   `scipy>=1.11.1`
-   `pandas>=2.0.3` (required for [Kratos solver wrapper](coupling_components/solver_wrappers/kratos_structure/kratos_structure.md))
-   `matplotlib>=3.7.2` (recommended)

We recommend Anaconda 2023.09 or newer. Older versions of the software might also work, but are not tested.


### Installation procedure

CoCoNuT does not need to be compiled, hence installation is straightforward. 
The source code can be downloaded as a zip file, or cloned directly from GitHub. For users that have no experience with Git or GitHub, we recommend the first option. The second option makes it easier to update the software and contribute to the code. 

*Option 1: download zip*

-   Download the source code from [GitHub](https://github.com/pyfsi/coconut).
-   Unzip to a folder *`coconut`*. If the folder is unzipped in Windows, some file permissions may change and some tests or examples may not run out of the box. 

*Option 2: clone source*

-   Choose or create a directory to install CoCoNuT. 
-   Move to this directory. 
-   Clone the GitHub repository with SSH or HTTPS by executing one of the following commands. 

    -   With SSH:
    
    ```
    git clone git@github.com:pyfsi/coconut.git
    ```
   
    -   With HTTPS:
         
    ```
    git clone https://github.com/pyfsi/coconut.git
    ```

After the code has been downloaded or cloned, the *`coconut`* folder must be added to the user's Python path. 
For example, with a folder structure like

```
/some/absolute/path/
    coconut/
        coupling_components/
        data_structure/
        ...
        README.md
```

*`coconut`* can be added to your Python path by executing the following line:

```bash
export PYTHONPATH=/some/absolute/path:$PYTHONPATH
```

This line can also be added to your *`.bashrc`* file.

### Checking the solver modules 

Before using CoCoNuT, it is necessary to adapt some system specific commands in the *`solver_modules.py`* file in the *`coconut`* folder.
This file has the commands to load solver modules in separate environments when running a case, to avoid conflicts. As these commands are system specific, it is important to check this file before testing CoCoNuT. 
The file contains a nested dictionary `solver_load_cmd_dict`, which has keys such as `ugent_cluster_SL6.3` or `ugent_cluster_CO7` denoting the machine on which CoCoNuT is installed.
In their turn, each of these dictionaries contains keys for all solvers that are available on that machine and can be used in CoCoNuT. 
The values are strings containing terminal commands to load the software, thus setting the environment which allows running the solver. 
For example, on the UGent cluster, the [Lmod system](https://lmod.readthedocs.io/en/latest/) is used, but there is no general guideline on how to make the solvers' software available as long as it is compatible with your system's command-line-interface. 
If multiple commands are needed, they should appropriately be separated within the string. For example in a Linux terminal the semicolon (;) or double ampersand (&&) can be used.
Since `machine_name` is set to `ugent_cluster_CO7`, this dictionary is used by default.
In case your system differs from the `ugent_cluster_CO7` settings, it is advised to add your own internal dictionary to `solver_load_cmd_dict` and provide this key to `machine_name`.
If a solver module is not present on your system the key should be removed. If a solver module is always present, i.e. no module load command or similar action is needed, an empty string should be given as value.
When CoCoNuT tries to use a solver module that is not present in the `solver_load_cmd_dict` or that has the wrong value, an error will be raised.
The standard output and error are redirected to a file named *`solver_load_cmd.log`*.

### Quick test

We recommend to run the unit tests at the end of the installation, to make sure that everything works. 

-   Ensure that *`coconut`* is included in your Python path.
-   Move to the *`coconut/tests`* directory. 
-   Run the _fast_ unit tests (excluding the non-Pythonic solvers) by executing the following line:
```bash
python3 run_tests.py
```
-   Or choose to include the solver wrapper tests by using the keyword `-all` as follows
```bash
python3 run_tests.py -all
```
More information on running tests can be found [here](tests/tests.md).


## Getting started

Once the CoCoNuT package has been successfully installed, it is time to run a first coupled simulation. For this purpose, we give a step-by-step guide of an example case included in the source code.

In this example the fluid-structure interaction (FSI) problem of a pressure wave propagating through an elastic tube in incompressible flow is calculated [[1](#1), [2](#2)]. For both the flow and structural solver, we use 1D Python-solvers that are included in CoCoNuT. This has the advantage that no external single-physics solvers must be installed for this example. Furthermore, the 1D solvers are very fast, so that a full transient FSI calculation can be done in this example. Other example cases in the source code solve the same FSI problem with ANSYS Fluent or OpenFOAM as flow solver and Abaqus or Kratos as structural solver. 

We start by creating a variable `COCO` in which we can store the path to the folder in which CoCoNuT is installed. We will use this variable to avoid any confusion about relative or absolute paths in this tutorial. Using the example installation location from above:

```bash
COCO=/some/absolute/path
```

We can now navigate to the folder of the example we will simulate. 
```bash
cd $COCO/coconut/examples/tube/tube_flow_tube_structure/
```
This folder serves as main directory to set up and run the FSI simulation from in CoCoNuT. The file *`parameters.json`* will be used to run the actual FSI simulation, but we will come back to that later. 
First we must set up both single-physics solvers separately. This setup is typically done outside CoCoNuT by the user, as it is solver and case specific. 
In this case we provide a script *`setup_case.py`* that sets up both solvers using the files in the folder *`../setup_files`*. When the script is run with

```bash
python3 setup_case.py
```

new folders *`CFD`* and *`CSM`* appear, as well as the file *`run_simulation.py`*. The *`CFD`* folder contains all files required to start a simulation of the flow in the tube. 
Analogously, the *`CSM`* folder contains all files required to start a simulation of the tube structure.

We can now start the FSI simulation in CoCoNuT by running the Python file *`run_simulation.py`*:

```bash
python3 run_simulation.py
```

The simulation should start, first printing the CoCoNuT ASCII-banner and some information about the settings of the FSI simulation. Then the simulation itself starts: in each time step, the residual is given for every coupling iteration. When the simulation has finished, a summary about the computational effort is printed.

Let us now take a closer look at the two files that are used to run CoCoNuT. 
The Python file *`run_simulation.py`* typically does not have to be adapted by the user. Its task is to read in the settings file *`parameters.json`* and launch a simulation using those settings. 
The file *`parameters.json`* is a collection of settings that is written in [JSON format](https://www.json.org/json-en.html). JSON is a language-independent text format that is easy to read and write, and is used for data-exchange. 
It consists mainly of key-value pairs, and can hence be easily converted to a (nested) Python dictionary. While the keys are always strings, the values can be strings, numbers, arrays, booleans or nested JSON objects (nested dictionaries).
Before you read on, it can be useful to familiarize yourself with the JSON syntax. In what follows, we will use Python terminology (dictionary, list, boolean, etc...) to refer to the structure and the values in the JSON file. 

The JSON file is built up in a hierarchical way that represents the objects created in the CoCoNuT simulation. At the highest level, the dictionary contains two keys: `settings` and `coupled_solver`. 
The value given to the `settings` key is a nested dictionary, which contains a single key-value pair that sets the number of time steps to be simulated. 
The value given to the `coupled_solver` key is a special dictionary, because it has the `type` key. CoCoNuT will generate an object of the specified type, namely `coupled_solvers.iqni`. This refers to the class defined in the file *`$COCO/coconut/coupling_components/coupled_solvers/iqni.py`*: the `CoupledSolverIQNI` class. 
Note that the value in `type` always refers to a file located in *`$COCO/coconut/coupling_components`*. 
The dictionary under `settings` is used to initialize an instance of this class. In this case the initial time `timestep_start`, the time step `delta_t` and some other parameters must be given. The coupled solver is the main class that determines how the two single-physics solvers are coupled. 
The dictionary that is given to the `coupled_solver` key contains next to `type` and `settings` three other key-value pairs. These will generate other objects: the fact that they are given in the `coupled_solver` dictionary means that these objects will be created by the coupled solver object.

`predictor` will generate an object of the `PredictorLinear` class found in the file *`$COCO/coconut/coupling_components/predictors/linear.py`*. This class requires no additional settings for its initialization. The predictor object is used to extrapolate the solution to the next time step. 

`convergence_criterion` will generate an object of the `ConvergenceCriterionOr` class found in the file *`$COCO/coconut/coupling_components/convergence_criteria/or.py`*, using the given `settings` for its initialization. The convergence criterion is used to determine when CoCoNuT should move to the next time step. In this case the *or* criterion is used, which signals convergence when one or both underlying criteria are satisfied. These underlying criteria are instances of the `ConvergenceCriterionIterationLimit` and `ConvergenceCriterionRelativeNorm` classes defined in respectively *`$COCO/coconut/coupling_components/convergence_criteria/iteration_limit.py`* and *`$COCO/coconut/coupling_components/convergence_criteria/relative_norm.py`*. 
This means that CoCoNuT will move to the next time step after 15 iterations or when the 2-norm of the residual has decreased six orders of magnitude.

`solver_wrappers` is a list of two solver wrapper objects, which will communicate with the two single-physics solvers, in this case the 1D flow solver and the 1D structural solver. The first dictionary in the list will generate an instance of the `SolverWrapperTubeFlow` class found in *`$COCO/coconut/coupling_components/solver_wrappers/python/tube_flow_solver.py`*. An important setting to generate this object is the `working_directory`, which refers to the folder *`CFD`* that we created with the case files of the flow solver. All files written by the flow solver will also appear in this folder.
We would now expect the second dictionary to generate a solver wrapper to communicate with the structural solver, i.e. an instance of the `SolverWrapperTubeStructure` class found in *`$COCO/coconut/coupling_components/solver_wrappers/python/tube_structure_solver.py`*. This is not the case however: the flow and structural solvers typically use a different geometrical discretization (computational grid or mesh), hence they cannot readily be coupled in CoCoNuT. To overcome this issue, we put a layer of mapping around one of the solver wrappers. This is done with the `SolverWrapperMapped` class found in *`$COCO/coconut/coupling_components/solver_wrappers/mapped.py`*. The *mapped* solver wrapper interpolates all data flowing between the coupled solver and the real solver wrapper. The mapped solver wrapper itself contains three objects: the actual solver wrapper (`SolverWrapperTubeStructure` class), and mappers for respectively the input and the output of the solver wrapper (both `MapperInterface` class, found in *`$COCO/coconut/coupling_components/mappers/interface.py`*). 

The concept of the mapped solver wrapper illustrates the modularity of CoCoNuT. As far as the coupled solver is concerned, the mapped solver wrapper acts exactly as a real solver wrapper. The real solver wrapper does not know about the mapping at all: it acts as if it directly communicates with the coupled solver. Furthermore, the interpolation method can be easily changed by swapping the mappers in the mapped solver wrapper: the current linear interpolation scheme can for example be replaced by a radial basis scheme by changing `mappers.linear` to `mappers.radial_basis`. 

Now try to change some settings in the JSON file, such as the mappers, the time step or the maximum number of coupling iterations, and rerun the coupled simulation.

After a simulation is finished, it can be useful to inspect or visualize the output quantities (i.e. displacement, pressure and in general also shear).
CoCoNuT has some built-in tools to do just that described in the [post-processing documentation](examples/post_processing/post_processing.md).
For the FSI-simulation we have just performed, an example script is present in *`$COCO/coconut/examples/post_processing/`*.
It requires the `write_results` setting in the `coupled_solver` part of the JSON-file to be set on a non-zero integer, which is for all examples done by default.
By running this example script *`animate_example.py`*, we will generate several animations:

```bash
python3 $COCO/coconut/examples/post_processing/animate_example.py
```

Animations of the displacement, pressure and coordinates (with equally scaled axes) will be shown.
The first two are shown below.

<p align="center">
  <img alt="Displacement animation" src="https://raw.githubusercontent.com/pyfsi/coconut/master/docs/images/displacement.gif" width="49%">
  <img alt="Pressure animation" src="https://raw.githubusercontent.com/pyfsi/coconut/master/docs/images/pressure.gif" width="49%">
</p>


## Overview of the code

The CoCoNuT package consists of 5 main folders: *`coupling_components`*, *`data_structure`*, *`docs`*, *`examples`* and *`tests`*. To give a general understanding of how the code is structured, we give a brief description of the purpose of each folder. The documentation website mirrors this folder structure and the folder names below link to the corresponding page.


### [*`coupling_components`*](coupling_components/coupling_components.md)

This folder contains the basic building blocks of CoCoNuT, which can be used to set up a coupled simulation. This includes among others the solver wrappers, to communicate with single-physics solvers, and the mappers, which provide interpolation between non-conforming meshes present in the different single-physics solvers.

### [*`data_structure`*](data_structure/data_structure.md)

This folder contains the data structure that is used internally in CoCoNuT to store and pass around information obtained from the single-physics solvers. The data structure relies on NumPy arrays for efficient storage and manipulation of data.


### [*`docs`*](docs/docs.md)

This folder serves to automatically generate the documentation website, based on the MarkDown documentation files that are present throughout the code. 


### [*`examples`*](examples/examples.md)

This folder contains examples of several fluid-structure interaction cases, which can serve as starting point for settings up the user's own simulation. They also provide insight into the capabilities of CoCoNuT.


### [*`tests`*](tests/tests.md)

This folder contains the unit tests. These are created for each piece of code that is added to CoCoNuT and are run regularly, to avoid bugs. 

## References
<a id="1">[1]</a> 
[Delaissé N., Demeester T., Haelterman R. and Degroote J., "Quasi-Newton methods for partitioned simulation of fluid-structure interaction reviewed in the generalized Broyden framework", Archives of Computational Methods in Engineering, vol. 30, pp. 3271-3300, 2023.](https://doi.org/10.1007/s11831-023-09907-y)

<a id="2">[2]</a> 
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Surrogate-based acceleration of quasi-Newton techniques for fluid-structure interaction simulations", Computers & Structures, vol. 260, pp. 106720, 2022.](http://hdl.handle.net/1854/LU-8728347)
