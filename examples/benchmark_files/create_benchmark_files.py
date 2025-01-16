import coconut
from coconut import tools

import json
import os
import shutil
import time
from shutil import copytree, copy
from os.path import join, isdir

"""
This file creates benchmark pickle files used by the script 'evaluate_examples.py'.

Running this file creates pickle files for the examples in the dictionary 'examples'.
The key is the folder name of the example and the corresponding value is an int indicating the number of 
time steps the example should be run.
In order to exclude some examples, comment out the corresponding key-value pair in the dictionary 'examples'.

To include another example, simply add a new key-value pair to the dictionary 'examples'.
"""

# set paths
execute_path = os.path.realpath(os.path.dirname(__file__))
tmp_path = join(execute_path, 'tmp')
examples_path = os.path.realpath(join(os.path.dirname(__file__), '..'))
benchmark = join(examples_path, 'benchmark_files')

# example and tuple with number of time steps and optional additional files
# comment out examples that you do not want to make benchmark files for
examples = {
    'breaking_dam/fluent2d_abaqus2d': 2,
    'breaking_dam/fluent2d_kratos_structure2d': 2,
    'lid_driven_cavity/fluent2d_kratos_structure2d': 2,
    'lid_driven_cavity/openfoam2d_kratos_structure2d': 2,
    'tube/fluent2d_abaqus2d': 2,
    'tube/fluent2d_abaqus2d_steady': 1,
    'tube/fluent2d_abaqus2d_surrogate': 2,
    'tube/fluent2d_tube_structure': 2,
    'tube/fluent3d_abaqus2d': 2,
    'tube/fluent3d_abaqus3d': 2,
    'tube/fluent3d_kratos_structure3d': 2,
    'tube/openfoam3d_abaqus3d': 2,
    'tube/openfoam3d_kratos_structure3d': 2,
    'tube/tube_flow_abaqus2d': 2,
    'tube/tube_flow_tube_ringmodel': 10,
    'tube/tube_flow_tube_structure': 10,
    'tube/tube_flow_tube_structure_analytical': 10,
    'tube/tube_flow_tube_structure_surrogate': 10,
    'turek/fluent2d_abaqus2d': 2,
    'turek/fluent2d_abaqus2d_steady': 1,
    'turek/openfoam2d_kratos_structure2d': 2,
}


def set_up():
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.mkdir(tmp_path)
    shutil.copy(join(examples_path, 'run_simulation.py'), tmp_path)


def clean_up():
    time.sleep(0.1)
    shutil.rmtree(tmp_path)


def create_benchmark(example, number_of_timesteps):
    # set paths
    tmp_example_path = join(tmp_path, example)

    # copy setup_files folder to tmp
    case = example.split('/')[0]
    if not isdir(join(tmp_path, case)):
        copytree(join(examples_path, case, "setup_files"), join(tmp_path, case, "setup_files"))

    # copy example folder to tmp
    os.makedirs(tmp_example_path)
    for file in ('parameters.json', 'setup_case.py'):
        copy(join(examples_path, example, file), tmp_example_path)

    # go to this example directory
    os.chdir(tmp_example_path)

    # perform set up
    tools.import_module('setup_case', join(tmp_example_path, 'setup_case.py'))

    # read parameters and limit number of time steps
    parameter_file_name = "parameters.json"
    with open(join(tmp_example_path, parameter_file_name), 'r') as parameter_file:
        parameters = json.load(parameter_file)
    parameters['settings']['number_of_timesteps'] = number_of_timesteps
    parameters['coupled_solver']['settings']['write_results'] = True
    parameters['coupled_solver']['settings']['case_name'] = 'case'
    parameters['coupled_solver']['settings']['anonymous'] = True

    # perform simulation
    simulation = coconut.Analysis(parameters)
    simulation.run()

    # return to initial execution path
    os.chdir(execute_path)

    # move and rename results data
    if not isdir(join(benchmark, case)):
        os.mkdir(join(benchmark, case))
    shutil.move(join(tmp_example_path, 'case_results.pickle'), join(benchmark, f'{example}.pickle'))


if __name__ == '__main__':
    set_up()
    for example, number_of_timesteps in examples.items():
        with tools.quick_timer(example):
            create_benchmark(example, number_of_timesteps)
    clean_up()
