from coconut import data_structure
from coconut.coupling_components.analysis import Analysis


if __name__ == '__main__':
    from sys import argv

    # Check number of command line arguments
    if len(argv) != 2:
        err_msg = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python run_simulation.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    # Import data structure
    parameter_file_name = argv[1]

    # Import parameters using the data structure
    with open(parameter_file_name, 'r') as parameter_file:
        parameters = data_structure.Parameters(parameter_file.read())

    simulation = Analysis(parameters)
    simulation.Run()
