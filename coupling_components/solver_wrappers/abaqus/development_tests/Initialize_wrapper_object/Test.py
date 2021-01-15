from coconut import data_structure
from coconut.tools import create_instance

from sys import argv


# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python Test.py <parameter-file>.json"\n'
    raise Exception(err_msg)

# Import data structure
parameter_file_name = argv[1]

# Import parameters using the data structure
with open(parameter_file_name, 'r') as parameter_file:
    parameters = data_structure.Parameters(parameter_file.read())

print("Going to create AbaqusSolver")
AbaqusSolver0  = create_instance(parameters["solver_wrappers"][0])
print("AbaqusSolver0 created")