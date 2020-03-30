# CoCoNuT installation

## For users

Requirements:

-   Python 3.6+

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




## For developers

Requirements:

-   Python 3.6+
-   Git

Installation:

-   Choose or create a directory to install CoCoNuT. 
-   Move to this directory. 
-   Execute the following git command to clone the repository.

```bash
git clone git@github.com:pyfsi/coconut.git
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






