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
This line can also be added to you `.bashrc` file.




## For developers

Requirements:

-   Python 3.6+
-   Git

Installation:

-   Choose or create a directory to install CoCoNuT.
-   Copy the code below to a bash-script and adapt the parameter `absolute_path` to the above directory.
-   Run bash-script.

```bash
#!/bin/bash

# add absolute path to parent directory of coconut
absolute_path=/some/absolute/path

# configure repository
cd $absolute_path
git clone git@github.com:pyfsi/coconut.git
cd coconut
git remote add upstream git@github.com:pyfsi/coconut.git

# add coconut to Python path --> put this line also in your .bashrc
export PYTHONPATH=$absolute_path:$PYTHONPATH

# run tests
cd $absolute_path/coconut/tests
python test_coconut.py

# run an FSI test example
cd $absolute_path/coconut/test_examples/tube_tube_flow_tube_structure
chmod 700 setup_tube_flow.sh setup_tube_structure.sh
./setup_tube_flow.sh
./setup_tube_structure.sh
python MainKratos.py project_parameters_mapped.json  # > $absolute_path/tube.log
```




