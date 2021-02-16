from coconut import data_structure
from coconut.tools import create_instance

import numpy as np
from sys import argv
import os
import json


def print_colored(string, color):
    if color=='green':
        print('\x1b[0;30;42m'+string+'\x1b[0m')
    elif color=='orange':
        print('\x1b[0;30;43m' + string + '\x1b[0m')
    elif color=='red':
        print('\x1b[0;30;41m' + string + '\x1b[0m')
    else:
        print(string+f'(color {color} not implemented)')


# Check number of command line arguments
if len(argv) != 2:
    err_msg = 'Wrong number of input arguments!\n'
    err_msg += 'Use this script in the following way:\n'
    err_msg += '    "python Test.py <parameter-file>.json"\n'
    raise Exception(err_msg)

# Import data structure
parameter_file_name = argv[1]

# Import parameters using the data structure
with open('parameters.json') as parameter_file:
    parameters = json.load(parameter_file)

# Create the solver (__init__)
print("Creating an AbaqusSolver")
AbaqusSolver0 = create_instance(parameters["solver_wrappers"][0])
print_colored("AbaqusSolver0 created", "green")

# Assign loads to the Input-Nodes
# give value to DISPLACEMENT variable
mp = AbaqusSolver0.model['BEAMINSIDEMOVING_load_points']  # interface input modelpart
pressure = vars(data_structure)['PRESSURE']
traction = vars(data_structure)['TRACTION']
p = 10000
for node in mp.Nodes:
    # Domain extends from Y -0.025 to 0.025, default x-position is 0.005
    node.SetSolutionStepValue(pressure, 0, p)
    node.SetSolutionStepValue(traction, 0, [0, 0, 0])
print(f"Assigned uniform pressure ({p} Pa) and 0 traction at the interface ")

AbaqusSolver0.initialize()

# Step 1, Coupling 1
AbaqusSolver0.initialize_solution_step()
AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())

os.system("cp -r CSM/CSM_Time1.odb CSM/CSM_Time1_Iter1.odb")

# Step 1, Coupling 2
AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())
AbaqusSolver0.finalize_solution_step()

#Step 2, Coupling 1
AbaqusSolver0.initialize_solution_step()
AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())
AbaqusSolver0.finalize_solution_step()

#Iterate until deformation is approximately steady
mp_out = AbaqusSolver0.model['BEAMINSIDEMOVING_nodes']  # interface input modelpart
displacement = vars(data_structure)['DISPLACEMENT']

n_out = mp_out.NumberOfNodes()
prev_displacement = np.zeros((n_out, 3))*0.
diff = np.zeros((n_out,3))*0.
tol = 1e-06
diffMax = 1000
while diffMax > tol:
    AbaqusSolver0.initialize_solution_step()
    AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())
    AbaqusSolver0.finalize_solution_step()
    diffMax = 0
    for node in mp_out.Nodes:
        diff = np.linalg.norm(np.array(node.GetSolutionStepValue(displacement))-prev_displacement[int(node.Id), :])
        prev_displacement[int(node.Id), :] = np.array(node.GetSolutionStepValue(displacement))
        if diff > diffMax:
            diffMax = diff
    print(diffMax)


AbaqusSolver0.finalize()

print_colored("Finished",'green')