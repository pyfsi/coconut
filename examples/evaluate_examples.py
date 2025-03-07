import coconut
from coconut import tools

import numpy as np
import pickle
import unittest
import json
import os
import shutil
from unittest import TestSuite, TestLoader
from os.path import join, isdir

"""
This file runs the example cases and compares the results and convergence history, with previously save benchmark
pickle files.

Run this file with the unittest framework (see unittest documentation for command-line flags):
    python3 evaluate_examples.py
Or just as a Python file.
In order to exclude some examples from running, comment out the corresponding class in the tuple 'test_cases'.

To include another example, create a new class analogous to the existing classes and add the class name to the tuple
'test_cases'. The class variable 'example' is the folder name and the class variable 'number_of_timesteps' refers to the 
number of time steps that the example should be run (the benchmark file should have run for at least this number of 
time steps).
There are two other class variables; 'additional_files' is a list of additional files in the example required for the 
simulation besides 'setup_case.py' and 'parameters.json', and 'compare_data' is a boolean which determines whether the
comparison with a benchmark pickle file should be made.
Once this class is made, you may need to create a benchmark file for this new example using the corresponding script.
"""

# set paths
execute_path = os.path.realpath(os.path.dirname(__file__))
tmp_path = join(execute_path, 'tmp')
examples_path = os.path.realpath(os.path.dirname(__file__))
bench_mark = join(examples_path, 'benchmark_files')


def setUpModule():
    if os.path.exists(tmp_path):
        shutil.rmtree(tmp_path)
    os.mkdir(tmp_path)
    shutil.copy(join(examples_path, 'run_simulation.py'), tmp_path)


def tearDownModule():
    tools.rm_timed(tmp_path)


class EvaluateExamples(unittest.TestCase):
    example = ''
    number_of_timesteps = 5
    additional_files = []
    compare_data = True
    atol_solution_x = 1e-11
    rtol_solution_x = 1e-4
    atol_solution_y = 1e-7
    rtol_solution_y = 1e-4
    atol_convergence = 1e-14
    rtol_convergence = 1e-4

    @classmethod
    def setUpClass(cls):
        # set paths
        tmp_example_path = join(tmp_path, cls.example)

        # copy setup_files folder to tmp
        case = cls.example.split('/')[0]
        if case == 'test_single_solver':
            case = 'tube'
        if not isdir(join(tmp_path, case)):
            shutil.copytree(join(examples_path, case, 'setup_files'), join(tmp_path, case, 'setup_files'))

        # copy example folder to tmp
        os.makedirs(tmp_example_path)
        for file in ['parameters.json', 'setup_case.py'] + cls.additional_files:
            shutil.copy(join(examples_path, cls.example, file), tmp_example_path)

        # go to this example directory
        os.chdir(tmp_example_path)

        # perform set up
        tools.import_module('setup_case', join(tmp_example_path, 'setup_case.py'))

        # read parameters and limit number of time steps
        parameter_file_name = 'parameters.json'
        with open(join(tmp_example_path, parameter_file_name), 'r') as parameter_file:
            parameters = json.load(parameter_file)
        parameters['settings']['number_of_timesteps'] = cls.number_of_timesteps
        parameters['coupled_solver']['settings']['write_results'] = True
        case_name = cls.example.replace('/', '_')
        parameters['coupled_solver']['settings']['case_name'] = case_name

        # perform simulation
        simulation = coconut.Analysis(parameters)
        simulation.run()

        # return to initial execution path
        os.chdir(execute_path)

        if cls.compare_data:
            # read results data
            with open(join(tmp_example_path, f'{case_name}_results.pickle'), 'rb') as f:
                cls.results = pickle.load(f)

            # read benchmark data
            with open(join(bench_mark, f'{cls.example}.pickle'), 'rb') as f:
                cls.benchmark = pickle.load(f)

    def test_solution(self):
        if self.compare_data:
            for interface, rtol, atol in (('solution_x', self.rtol_solution_x, self.atol_solution_x),
                                          ('solution_y', self.rtol_solution_y, self.atol_solution_y)):
                np.testing.assert_allclose(self.benchmark[interface][:, :self.number_of_timesteps + 1],
                                           self.results[interface][:, :self.number_of_timesteps + 1],
                                           rtol=rtol, atol=atol)
        else:
            self.skipTest(f'Solution not compared for {self.__class__.__name__}')

    def test_convergence_history(self):
        if self.compare_data:
            for time_step in range(self.number_of_timesteps):
                np.testing.assert_allclose(np.array(self.benchmark['residual'][time_step]),
                                           np.array(self.results['residual'][time_step]),
                                           rtol=self.rtol_convergence, atol=self.atol_convergence)
        else:
            self.skipTest(f'Convergence history not compared for {self.__class__.__name__}')


class TestBreakingDamFluent2DAbaqus2D(EvaluateExamples):
    example = 'breaking_dam/fluent2d_abaqus2d'
    number_of_timesteps = 2


class TestBreakingDamFluent2DKratosStructure2D(EvaluateExamples):
    example = 'breaking_dam/fluent2d_kratos_structure2d'
    number_of_timesteps = 2


class TestLidDrivenCavityFluent2DKratosStructure2D(EvaluateExamples):
    example = 'lid_driven_cavity/fluent2d_kratos_structure2d'
    number_of_timesteps = 2


class TestLidDrivenCavityOpenFOAM2DKratosStructure2D(EvaluateExamples):
    example = 'lid_driven_cavity/openfoam2d_kratos_structure2d'
    number_of_timesteps = 2


class TestTestSingleSolver(EvaluateExamples):
    example = 'test_single_solver'
    number_of_timesteps = 2
    additional_files = ['dummy_solver.py']
    compare_data = False


class TestTubeFluent2DAbaqus2D(EvaluateExamples):
    example = 'tube/fluent2d_abaqus2d'
    number_of_timesteps = 2


class TestTubeFluent2DAbaqus2DSteady(EvaluateExamples):
    example = 'tube/fluent2d_abaqus2d_steady'
    number_of_timesteps = 1


class TestTubeFluent2DAbaqus2DSurrogate(EvaluateExamples):
    example = 'tube/fluent2d_abaqus2d_surrogate'
    number_of_timesteps = 2


class TestTubeFluent2DTubeStructure(EvaluateExamples):
    example = 'tube/fluent2d_tube_structure'
    number_of_timesteps = 2


class TestTubeFluent3DAbaqus2D(EvaluateExamples):
    example = 'tube/fluent3d_abaqus2d'
    number_of_timesteps = 2
    atol_solution_y = 1e-6


class TestTubeFluent3DAbaqus3D(EvaluateExamples):
    example = 'tube/fluent3d_abaqus3d'
    number_of_timesteps = 2


class TestTubeFluent3DKratosStructure3D(EvaluateExamples):
    example = 'tube/fluent3d_kratos_structure3d'
    number_of_timesteps = 2


class TestTubeOpenFOAM3DAbaqus3D(EvaluateExamples):
    example = 'tube/openfoam3d_abaqus3d'
    number_of_timesteps = 2
    atol_solution_y = 3e-4


class TestTubeOpenFOAM3DKratosStructure3D(EvaluateExamples):
    example = 'tube/openfoam3d_kratos_structure3d'
    number_of_timesteps = 2
    atol_solution_y = 3e-4


class TestTubeTubeFlowAbaqus2D(EvaluateExamples):
    example = 'tube/tube_flow_abaqus2d'
    number_of_timesteps = 2


class TestTubeTubeFlowTubeRingmodel(EvaluateExamples):
    example = 'tube/tube_flow_tube_ringmodel'
    number_of_timesteps = 5
    atol_solution_x = 1e-12
    rtol_solution_x = 0
    atol_solution_y = 1e-5
    rtol_solution_y = 0
    atol_convergence = 1e-9
    rtol_convergence = 0


class TestTubeTubeFlowTubeStructure(EvaluateExamples):
    example = 'tube/tube_flow_tube_structure'
    number_of_timesteps = 5
    atol_solution_x = 1e-12
    rtol_solution_x = 0
    atol_solution_y = 1e-5
    rtol_solution_y = 0
    atol_convergence = 1e-9
    rtol_convergence = 0


class TestTubeTubeFlowTubeStructureAnalytical(EvaluateExamples):
    example = 'tube/tube_flow_tube_structure_analytical'
    number_of_timesteps = 5
    atol_solution_x = 1e-12
    rtol_solution_x = 0
    atol_solution_y = 1e-5
    rtol_solution_y = 0
    atol_convergence = 1e-9
    rtol_convergence = 0


class TestTubeTubeFlowTubeStructureSurrogate(EvaluateExamples):
    example = 'tube/tube_flow_tube_structure_surrogate'
    number_of_timesteps = 5
    atol_solution_x = 1e-12
    rtol_solution_x = 0
    atol_solution_y = 1e-5
    rtol_solution_y = 0
    atol_convergence = 1e-9
    rtol_convergence = 0


class TestTurekFluent2DAbaqus2D(EvaluateExamples):
    example = 'turek/fluent2d_abaqus2d'
    number_of_timesteps = 2


class TestTurekFluent2DAbaqus2DSteady(EvaluateExamples):
    example = 'turek/fluent2d_abaqus2d_steady'
    number_of_timesteps = 1


class TestTurekOpenFOAM2DKratosStructure2D(EvaluateExamples):
    example = 'turek/openfoam2d_kratos_structure2d'
    number_of_timesteps = 2


# comment out examples that you do not want to evaluate
test_cases = (
    TestBreakingDamFluent2DAbaqus2D,
    TestBreakingDamFluent2DKratosStructure2D,
    TestLidDrivenCavityFluent2DKratosStructure2D,
    TestLidDrivenCavityOpenFOAM2DKratosStructure2D,
    TestTestSingleSolver,
    TestTubeFluent2DAbaqus2D,
    TestTubeFluent2DAbaqus2DSteady,
    TestTubeFluent2DAbaqus2DSurrogate,
    TestTubeFluent2DTubeStructure,
    TestTubeFluent3DAbaqus2D,
    TestTubeFluent3DAbaqus3D,
    TestTubeFluent3DKratosStructure3D,
    TestTubeOpenFOAM3DAbaqus3D,
    TestTubeOpenFOAM3DKratosStructure3D,
    TestTubeTubeFlowAbaqus2D,
    TestTubeTubeFlowTubeRingmodel,
    TestTubeTubeFlowTubeStructure,
    TestTubeTubeFlowTubeStructureAnalytical,
    TestTubeTubeFlowTubeStructureSurrogate,
    TestTurekFluent2DAbaqus2D,
    TestTurekFluent2DAbaqus2DSteady,
    TestTurekOpenFOAM2DKratosStructure2D,
)


def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


def suite():
    suite = TestSuite()
    loader = TestLoader()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
