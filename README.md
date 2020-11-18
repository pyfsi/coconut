# CoCoNuT - Coupling Code for Numerical Tools


CoCoNuT is a Python package for partitioned multi-physics simulations, with a focus on fluid-structure interaction. 
Thanks to its fully modular approach, the package is versatile and easy to extend. It is available under the GPL-3.0 license. 


## Introduction

The *Coupling Code for Numerical Tools*&mdash;*CoCoNuT* in short&mdash;follows a partitioned approach to solve multi-physics problem: existing single-physics solvers are coupled through a Python interface. 
This has the advantage that dedicated, highly-optimized single-physics solvers can be used. 
To use the code with a new solver (open source or commercial), a so-called *solver-wrapper* is written to take care of the communication with CoCoNuT. 
All other CoCoNuT components, such as mapping for non-conformal meshes, are solver-independent and are easily swapped thanks to CoCoNuT's modular architecture. 

CoCoNuT is under active development by the [Fluid Mechanics research team](https://www.ugent.be/ea/eemmecs/en/research/stfes/flow) at Ghent University. 
Our specialization is partitioned fluid-structure interaction. 
We develop high-performance quasi-Newton algorithms to accelerate convergence of the coupled problem, and apply these techniques to diverse fluid-structure interaction applications such as wind turbines, tube bundles and flexible aircraft. 

The full documentation of the CoCoNuT package can be found at the [documentation website](https://pyfsi.github.io/coconut/).


## Installation

These instructions describe the setup of CoCoNuT on Linux. The package has not been tested on Windows or macOS, so compatibility is not guaranteed, although we do not expect major issues. 


### Requirements

-   Python 3.6+ 
-   NumPy and SciPy packages (we recommend Anaconda 2019.03+)


### Installation

CoCoNuT does not need to be compiled, hence installation is straightforward. 
The source code can be downloaded as a zip file, or cloned directly from GitHub. For users that have no experience with Git or GitHub, we recommend the first option. The second option makes it easier to update the software and contribute to the code. 

*Option 1: download zip*

-   Download the source code from [GitHub](https://github.com/pyfsi/coconut).
-   Unzip to a folder `coconut`. If the folder is unzipped in Windows, some of the file permissions may change and some tests or examples may not run out of the box. 

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

After the code has been downloaded or cloned, the `coconut` package must be added to the user's Python path. 
For example, with a folder structure like

```
/some/absolute/path/
    coconut/
        coupling_components/
        data_structure/
        ...
        README.md
```

`coconut` can be added to your Python path by executing the following line:

```bash
export PYTHONPATH=/some/absolute/path:$PYTHONPATH
```

This line can also be added to your `.bashrc` file. 


### Quick test

We recommend to run the unit tests at the end of the installation, to make sure that everything works. 

-   Ensure that `coconut` is included in your Python path.
-   Move to the `coconut` directory. 
-   Run the unit tests by executing the following line:

    ```bash
    sh run_coconut_tests.sh
    ```



## Getting started

> TODO: this section should contain a kind of tutorial that goes over one of the test cases and explains more or less step-by-step how you run it, and what happens when you run it.

> TODO: change 'test_examples' to 'examples' where necessary


Once the CoCoNuT package has been successfully installed, it is time to run a first example case. For this purpose, a step-by-step guide to run an example included in the source code is given.

In this example the fluid-structure interaction (FSI) problem of a pressure wave propagating through an elastic tube in incompressible flow is calculated. For both the flow and structure solver, we use 1D Python-solvers that are included in CoCoNuT. This has the advantage that no external single-physics solvers must be installed for this example. Furthermore, the 1D solvers are very fast, so that a full transient FSI calculation can be done in this example. 

We start by creating an environment variable `COCO` in which we can store the path to the folder in which CoCoNuT is installed. We will use this variable to avoid any confusion about relative or absolute paths in this tutorial. Using the example installation location from above:

```bash
COCO=/some/absolute/path
```

> TODO: perhaps this is overkill, as we don't use path very often actually

We can now navigate to the folder of the example we will simulate. 

```bash
cd $COCO/coconut/test_examples/tube_tube_flow_tube_structure/
```

This folder contains all the files required to set up and run the FSI simulation in CoCoNuT. The files `run_simulation.py` and `project_paramaters_mapped.json` will be used to run the actual FSI simulation, but we will come back to those later. 

First we must set up both single-physics solvers separately, starting with the flow solver. This setup is typically done outside of CoCoNuT by the user, as it is solver-specific, but in this case we provide a bash script `setup_tube_flow.sh` that uses the files in the folder `setup_tube_flow` to generate the case. When the script is run with

```bash
./setup_tube_flow.sh
```

a new folder `CFD` appears. This folder contains all files required to start a simulation of the flow in the tube. Analogously, we run 

```bash
./setup_tube_structure.sh
```

to generate the `CSM` folder, which contains all files required to start a simulation of the tube structure. 

We can now start the FSI simulation in CoCoNuT by running the Python file `run_simulation.py` with the settings file `project_paramaters_mapped.json` as argument:

```bash
python run_simulation.py project_parameters_mapped.json
```

The simulation should start, first printing the CoCoNuT ASCII-banner, then some information about the settings of the FSI simulation and finally the residuals of the coupling iterations per time step. 


> TODO:
> - run_simulation.py: this is standard file, typically not changed 
> - project_parameters.json: go over file, explain a whole bunch of settings; also suggest to change some parameters and rerun simulation; 
> - finally: show an exmple of post-processing? we'll see... 






## Overview of the code

The CoCoNuT package consists of 5 main folders: `coupling_components`, `data_structure`, `docs`, `examples` and `tests`. To give a general understanding of how the code is structured, we give a brief description of the purpose of each folder. The documentation website mirrors this folder structure.


### `coupling_components`

This folder contains the basic building blocks of CoCoNuT, which can be used to set up a coupled simulation. This includes among others the solver-wrappers, to communicate with single-physics solvers, and the mappers, which provide interpolation between non-conforming meshes present in the different single-physics solvers.

### `data_structure`

This folder contains the data structure that is used internally in CoCoNuT to store and pass around information obtained from the single-physics solvers. The data structure relies on NumPy arrays for efficient storage and manipulation of data.


### `docs`

This folder serves to automatically generate the documentation website, based on the MarkDown documentation files that are present throughout the code. 


### `examples`

This folder contains examples of several fluid-structure interaction cases, which can serve as starting point for settings up the user's own simulation. They also provide insight into the capabilities of CoCoNuT.


### `tests`

This folder contains the unit tests. These are created for each piece of code that is added to CoCoNuT and are run regularly, to avoid bugs. 
