# CoCoNuT installation

## For users

Requirements:

-   Python 3.6+ (Anaconda 2019.03+ recommended)

Installation:

-   Download the source code from [GitHub](https://github.com/pyfsi/coconut)
-   Unzip to a folder `coconut`.
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

Documentation:

All documentation is grouped on the [CoCoNut website][https://pyfsi.github.io/coconut/].


## For developers

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

Documentation:

All documentation is grouped on the [CoCoNut website][https://pyfsi.github.io/coconut/].






[https://pyfsi.github.io/coconut/]: https://pyfsi.github.io/coconut/