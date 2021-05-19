import coconut
from coconut import tools

import json
import os
import shutil
import subprocess
from shutil import copytree, copy
from os.path import join

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
    'tube_fluent2d_abaqus2d': 2,
    'tube_fluent2d_abaqus2d_steady': 1,
    'tube_fluent2d_tube_structure': 2,
    'tube_fluent3d_abaqus2d': 2,
    'tube_fluent3d_abaqus3d': 2,
    'tube_fluent3d_kratos_structure3d': 2,
    'tube_openfoam3d_abaqus3d': 2,
    'tube_openfoam3d_kratos_structure3d': 2,
    'tube_tube_flow_abaqus2d': 2,
    'tube_tube_flow_tube_ringmodel': 100,
    'tube_tube_flow_tube_structure': 100,
    'turek_fluent2d_abaqus2d': 2,
    'turek_fluent2d_abaqus2d_steady': 1
}


def set_up():
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.mkdir(tmp_path)
    copytree(join(examples_path, "setup_files"), join(tmp_path, "setup_files"))


def clean_up():
    shutil.rmtree(tmp_path)


def create_benchmark(example, number_of_timesteps):
    # set paths
    tmp_example_path = join(tmp_path, example)

    # copy example folder to tmp
    os.mkdir(tmp_example_path)
    for file in ('parameters.json', 'setup.py'):
        copy(join(examples_path, example, file), tmp_example_path)

    # go to this example directory
    os.chdir(tmp_example_path)

    # perform set up
    tools.import_module('setup', join(tmp_example_path, 'setup.py'))

    # read parameters and limit number of time steps
    parameter_file_name = "parameters.json"
    with open(join(tmp_example_path, parameter_file_name), 'r') as parameter_file:
        parameters = json.load(parameter_file)
    parameters['settings']['number_of_timesteps'] = number_of_timesteps
    parameters['coupled_solver']['settings']['save_results'] = True
    parameters['coupled_solver']['settings']['case_name'] = 'case'

    # perform simulation
    simulation = coconut.Analysis(parameters)
    simulation.run()

    # return to initial execution path
    os.chdir(execute_path)

    # move and rename results data
    shutil.move(join(tmp_example_path, 'case_results.pickle'), join(benchmark, f'{example}.pickle'))


if __name__ == '__main__':
    set_up()
    for example, number_of_timesteps in examples.items():
        with tools.quick_timer(example):
            create_benchmark(example, number_of_timesteps)
    clean_up()
