from coconut import data_structure
from coconut.coupling_components.tools import create_instance

import numpy as np
import json


with open('parameters.json') as parameter_file:
    parameters = json.load(parameter_file)

solver = create_instance(parameters)
settings = parameters['settings']

# steady test  *** not converted
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
                dy = (1 - np.cos(2 * np.pi * node.X0)) * 0.5 * 0.01  # this used node.X before
                node.SetSolutionStepValue(vars(data_structure)['DISPLACEMENT'], 0, [0., dy, 0.])

    solver.finalize_solution_step()
    solver.finalize()

# unsteady test  *** being converted
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
                    dy = (1 - np.cos(2 * np.pi * (node.X0 - timestep / 4 - iteration / 16))) * 0.5 * f  # this used node.X before
                    node.SetSolutionStepValue(vars(data_structure)['DISPLACEMENT'], 0, [0., dy, 0.])
        solver.finalize_solution_step()

    solver.finalize()