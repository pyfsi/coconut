import coconut
from coconut import tools

import json
import os
import shutil
import subprocess
from shutil import copy2
from os.path import join

# set paths
execute_path = os.path.realpath(os.path.dirname(__file__))
tmp_path = join(execute_path, 'tmp')
examples_path = os.path.realpath(join('..', os.path.dirname(__file__)))
bench_mark = join(examples_path, 'benchmark_files')


def set_up():
    os.mkdir(tmp_path)
    cmd = f'cp -r {join(examples_path, "setup_files")} {join(tmp_path, "setup_files")}'
    p = subprocess.Popen(cmd, executable='/bin/bash', cwd=tmp_path, shell=True)
    p.wait()


def clean_up():
    shutil.rmtree(tmp_path)


def create_benchmark(example, number_of_timesteps):
    # set paths
    tmp_example_path = join(tmp_path, example)

    # copy example folder to tmp
    os.mkdir(tmp_example_path)
    for file in ('parameters.json', 'setup.sh'):
        copy2(join(examples_path, example, file), tmp_example_path)

    # perform set up
    p = subprocess.Popen(join(tmp_example_path, 'setup.sh'), cwd=tmp_example_path, shell=True)
    p.wait()

    # read parameters and limit number of time steps
    parameter_file_name = "parameters.json"
    with open(join(tmp_example_path, parameter_file_name), 'r') as parameter_file:
        parameters = json.load(parameter_file)
    parameters['settings']['number_of_timesteps'] = number_of_timesteps

    # perform simulation
    os.chdir(tmp_example_path)
    simulation = coconut.Analysis(parameters)
    simulation.run()
    os.chdir(execute_path)

    # move and rename results data
    shutil.move(join(tmp_example_path, 'results.pickle'), join(bench_mark, f'{example}.pickle'))


# example and number of time steps
# comment out examples that you do not want to make benchmark files for
examples = {
    'test_single_solver': 2,
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
    'tube_tube_flow_tube_structure': 100
    }

set_up()
for example, number_of_timesteps in examples.items():
    with tools.quick_timer(example):
        create_benchmark(example, number_of_timesteps)
clean_up()
