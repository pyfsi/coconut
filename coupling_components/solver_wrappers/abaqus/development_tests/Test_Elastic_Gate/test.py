from coconut import data_structure
from coconut.coupling_components.tools import create_instance

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
with open('abaqus_params.json') as parameter_file:
    parameters = json.load(parameter_file)

# Create the solver (__init__)
print("Creating an AbaqusSolver")
AbaqusSolver0 = create_instance(parameters["solver_wrappers"][0])
print_colored("AbaqusSolver0 created", "green")

# Assign loads to the Input-Nodes
# give value to PRESSURE variable

interface_input = AbaqusSolver0.get_interface_input()
model_part_0_name = "BEAMINSIDEMOVING0_load_points"
model_part_1_name = "BEAMINSIDEMOVING1_load_points"
model_part_2_name = "BEAMINSIDEMOVING2_load_points"
model_part_0 = interface_input.get_model_part(model_part_0_name)
model_part_1 = interface_input.get_model_part(model_part_1_name)
model_part_2 = interface_input.get_model_part(model_part_2_name)

p = 10
t = np.array((0, 0, 0))

pressure = interface_input.get_variable_data(model_part_2_name, 'pressure')
pressure[:] = p
interface_input.set_variable_data(model_part_2_name, 'pressure', pressure)
traction = interface_input.get_variable_data(model_part_2_name, 'traction')
traction[:, :] = t
interface_input.set_variable_data(model_part_2_name, 'traction', traction)

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
interface_output = AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())
AbaqusSolver0.finalize_solution_step()

displacement = interface_output.get_variable_data('BEAMINSIDEMOVING2_nodes', "displacement")
mp_out = interface_output.get_model_part('BEAMINSIDEMOVING2_nodes')

tol = 1e-07

prev_displacement = displacement

max_diff = 1000
while max_diff > tol:
    AbaqusSolver0.initialize_solution_step()
    interface_output = AbaqusSolver0.solve_solution_step(AbaqusSolver0.get_interface_input())
    AbaqusSolver0.finalize_solution_step()
    displacement = interface_output.get_variable_data('BEAMINSIDEMOVING2_nodes', "displacement")
    diff = np.linalg.norm(displacement - prev_displacement, axis=1)
    max_disp = np.max(np.linalg.norm(displacement, axis=1))
    max_diff = np.max(diff)
    print(f'max_disp: {max_disp}')
    print(f'max_diff: {max_diff}')
    prev_displacement = displacement

AbaqusSolver0.finalize()

print_colored("Finished", 'green')