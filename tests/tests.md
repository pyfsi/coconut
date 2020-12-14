# Tests

## Running unittests
The unittests in CoCoNuT uses unittest module available in python. To run all the tests you have to first go to `tests/`
directory and run the following command in the terminal:

````
python -m unittest discover -bv
```` 
The discover keyword finds all the unittests writen in the files named `test_<some name>`. It will recursively find all the test files with this pattern for all the folders containing `__init__.py` file. Further  documentation on unittest discover can be found [here](https://docs.python.org/3/library/unittest.html).
By default, all the unitests except the ones for the solver wrappers will be run.


## Including unittests for solver-wrappers 
Solver wrapper provides an interface to an external solver that is coupled while performing coupled simulation. Since this depends on the availability of the solver in a system, the unittests for the solver wrappers are made optional. Including the unitests for solver wrappers in the list of unittests that need to be run involves following steps: 

-   Load the corresponding solvers.
-   Import the unittest classes for the solver wrappers in the python file, `tests/solver_wrappers/__init__.py`.
-   Add the unittest classes for the solver wrappers in the list called `tests_cases` in the python file, `tests/solver_wrappers/__init__.py`.

After following the above steps, if you now run all the tests using the command given in the previous section from the directory `tests/`, unitests for the added solver wrappers will be run along with all the default ones.
