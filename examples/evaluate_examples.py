import coconut

import numpy as np
import pickle
import unittest
import json
import os
import shutil
import subprocess
from unittest import TestSuite
from os.path import join


# set paths
execute_path = os.path.realpath(os.path.dirname(__file__))
tmp_path = join(execute_path, 'tmp')
examples_path = os.path.realpath(os.path.dirname(__file__))
bench_mark = join(examples_path, 'benchmark_files')


def setUpModule():
    os.mkdir(tmp_path)
    cmd = f'cp -r {join(examples_path, "setup_files")} {join(tmp_path, "setup_files")}'
    p = subprocess.Popen(cmd, executable='/bin/bash', cwd=tmp_path, shell=True)
    p.wait()


def tearDownModule():
    shutil.rmtree(tmp_path)


class EvaluateExamples(unittest.TestCase):
    example = ''
    number_of_timesteps = 5

    @classmethod
    def setUpClass(cls):
        # set paths
        tmp_example_path = join(tmp_path, cls.example)

        # copy example folder to tmp
        os.mkdir(tmp_example_path)
        for file in ('parameters.json', 'setup.sh'):
            shutil.copy2(join(examples_path, cls.example, file), tmp_example_path)

        # perform set up
        p = subprocess.Popen(join(tmp_example_path, 'setup.sh'), cwd=tmp_example_path, shell=True)
        p.wait()

        # read parameters and limit number of time steps
        parameter_file_name = "parameters.json"
        with open(join(tmp_example_path, parameter_file_name), 'r') as parameter_file:
            parameters = json.load(parameter_file)
        parameters['settings']['number_of_timesteps'] = cls.number_of_timesteps

        # perform simulation
        os.chdir(tmp_example_path)
        simulation = coconut.Analysis(parameters)
        simulation.run()
        os.chdir(execute_path)

        # read results data
        with open(join(tmp_example_path, 'results.pickle'), 'rb') as f:
            cls.results = pickle.load(f)

        # read benchmark data
        with open(join(bench_mark, f'{cls.example}.pickle'), 'rb') as f:
            cls.benchmark = pickle.load(f)

    def test_solution(self):
        for interface in ('solution_x', 'solution_y'):
            np.testing.assert_almost_equal(self.benchmark[interface][:, :self.number_of_timesteps],
                                           self.results[interface][:, :self.number_of_timesteps])

    def test_convergence_history(self):
        for time_step in range(self.number_of_timesteps):
            np.testing.assert_allclose(np.array(self.benchmark['residual'][time_step]),
                                       np.array(self.results['residual'][time_step]))


class TestTestSingleSolver(EvaluateExamples):
    example = 'test_single_solver'
    number_of_timesteps = 2


class TestTubeFluent2DAbaqus2D(EvaluateExamples):
    example = 'tube_fluent2d_abaqus2d'
    number_of_timesteps = 2


class TestTubeFluent2DAbaqus2DSteady(EvaluateExamples):
    example = 'tube_fluent2d_abaqus2d_steady'
    number_of_timesteps = 1


class TestTubeFluent2DTubeStructure(EvaluateExamples):
    example = 'tube_fluent2d_tube_structure'
    number_of_timesteps = 2


class TestTubeFluent3DAbaqus2D(EvaluateExamples):
    example = 'tube_fluent3d_abaqus2d'
    number_of_timesteps = 2


class TestTubeFluent3DAbaqus3D(EvaluateExamples):
    example = 'tube_fluent3d_abaqus3d'
    number_of_timesteps = 2


class TestTubeFluent3DKratosStructure3D(EvaluateExamples):
    example = 'tube_fluent3d_kratos_structure3d'
    number_of_timesteps = 2


class TestTubeOpenFOAM3DAbaqus3D(EvaluateExamples):
    example = 'tube_openfoam3d_abaqus3d'
    number_of_timesteps = 2


class TestTubeOpenFOAM3DKratosStructure3D(EvaluateExamples):
    example = 'tube_openfoam3d_kratos_structure3d'
    number_of_timesteps = 2


class TestTubeTubeFlowAbaqus2D(EvaluateExamples):
    example = 'tube_tube_flow_abaqus2d'
    number_of_timesteps = 2


class TestTubeTubeFlowTubeRingmodel(EvaluateExamples):
    example = 'tube_tube_flow_tube_ringmodel'
    number_of_timesteps = 100


class TestTubeTubeFlowTubeStructure(EvaluateExamples):
    example = 'tube_tube_flow_tube_structure'
    number_of_timesteps = 100


# comment out examples that you do not want to evaluate
test_cases = (
    TestTestSingleSolver,
    TestTubeFluent2DAbaqus2D,
    TestTubeFluent2DAbaqus2DSteady,
    TestTubeFluent2DTubeStructure,
    TestTubeFluent3DAbaqus2D,
    TestTubeFluent3DAbaqus3D,
    TestTubeFluent3DKratosStructure3D,
    TestTubeOpenFOAM3DAbaqus3D,
    TestTubeOpenFOAM3DKratosStructure3D,
    TestTubeTubeFlowAbaqus2D,
    TestTubeTubeFlowTubeRingmodel,
    TestTubeTubeFlowTubeStructure
)


def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite
