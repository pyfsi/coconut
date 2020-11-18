# CoCoNuT - Coupling Code for Numerical Tools


## Introduction

> TODO: this section should contain a general overview of how the coconut code works, without going into specifics, mostly the philosophy but keep it quite short.

## CoCoNuT installation

### For users

These instructions describe the setup of CoCoNuT on Linux. It has not been tested on Windows or macOS, so compatibility is not guaranteed, although it probably would work without major issues. 

Requirements:

-   Python 3.6+ (Anaconda 2019.03+ recommended)

Installation:

-   Download the source code from [GitHub](https://github.com/pyfsi/coconut)
-   Unzip to a folder `coconut`. If the folder is unzipped in Windows, some of the file permissions may change and some tests or test examples may not run out of the box. 
-   Add the parent folder of `coconut` to your Python path.

For example, with a folder structure like
```
/some/absolute/path/
    coconut/
        coupling_components/
        data_structure/
        ...
        README.md
```
`coconut` can be added to your Python path with:
```bash
export PYTHONPATH=/some/absolute/path:$PYTHONPATH
```
This line can also be added to your `.bashrc` file.

For information about the CoCoNuT package, we refer to the [documentation website][https://pyfsi.github.io/coconut/].


### For developers

Requirements:

-   Python 3.6+ (Anaconda 2019.03+ recommended)
-   Git

Installation:

-   Choose or create a directory to install CoCoNuT. 
-   Move to this directory. 
-   Clone the GitHub repository with SSH or HTTPS by executing one of the following commands. The choice depends on your git configuration.

    *   With SSH:
    
    ```
    git clone git@github.com:pyfsi/coconut.git
    ```
   
    *   With HTTPS
         
    ```
    git clone https://github.com/pyfsi/coconut.git
    ```
    
-   Move to the `coconut` directory. 
-   Load Anaconda module.
-   Run `run_coconut_tests.sh` to see if the installation works.

```bash
sh run_coconut_tests.sh
```
-   Add the parent folder of `coconut` to your Python path.

For example, with a folder structure like

```
/some/absolute/path/
    coconut/
        coupling_components/
        data_structure/
        ...
        README.md
```

`coconut` can be added to your Python path with:

```bash
export PYTHONPATH=/some/absolute/path:$PYTHONPATH
```

This line can also be added to your `.bashrc` file.

For information about the CoCoNuT package, we refer to the [documentation website][https://pyfsi.github.io/coconut/].

[https://pyfsi.github.io/coconut/]: https://pyfsi.github.io/coconut/


## Getting started

> TODO: this section should contain a kind of tutorial that goes over one of the test cases and explains more or less step-by-step how you run it, and what happens when you run it.



## Overview of the code

> TODO: this section gives an overview of the 5 main folders in coconut, the rest of the website follows the structure of these folders. Some details are given about what is in each folder (coupling_components, data_structure, docs, test_examples, tests).