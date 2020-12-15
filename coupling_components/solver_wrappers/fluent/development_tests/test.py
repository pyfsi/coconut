from coconut import data_structure
from coconut.coupling_components.tools import CreateInstance

import numpy as np
from sys import argv


# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python co_simulation_analysis.py <cosim-parameter-file>.json"\n'
    raise Exception(err_msg)


# Import data structure
parameter_file_name = argv[1]

# Import parameters using the data structure
with open(parameter_file_name, 'r') as parameter_file:
    parameters = data_structure.Parameters(parameter_file.read())

solver = CreateInstance(parameters['solver_wrappers'][0])


settings = parameters['solver_wrappers'][0]['settings']

# steady test
if 0:
    solver.initialize()
    solver.initialize_solution_step()

    interface_input = solver.get_interface_input()
    for iteration in range(3):
        iteration += 1
        print(f'\niteration {iteration}')
        solver.solve_solution_step(interface_input)
        interface_input = solver.get_interface_input()
        for key in settings['interface_input'].keys():
            for node in interface_input.model[key].Nodes:
                dy = (1 - np.cos(2 * np.pi * node.X)) * 0.5 * 0.01
                node.SetSolutionStepValue(vars(data_structure)['DISPLACEMENT'], 0, [0., dy, 0.])

    solver.finalize_solution_step()
    solver.finalize()

# unsteady test
else:
    solver.initialize()

    interface_input = solver.get_interface_input()
    for timestep in range(1, 5):
        f = 0.005 * (-1) ** (timestep + 1)
        f = 0.05
        solver.initialize_solution_step()
        for iteration in range(1, 3):
            solver.solve_solution_step(interface_input)
            interface_input = solver.get_interface_input()
            for key in settings['interface_input'].keys():
                for node in interface_input.model[key].Nodes:
                    dy = (1 - np.cos(2 * np.pi * (node.X - timestep / 4 - iteration / 16))) * 0.5 * f
                    node.SetSolutionStepValue(vars(data_structure)['DISPLACEMENT'], 0, [0., dy, 0.])
        solver.finalize_solution_step()

    solver.finalize()