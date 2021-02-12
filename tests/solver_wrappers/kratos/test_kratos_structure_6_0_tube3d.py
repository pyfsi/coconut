from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.interface import Interface

import numpy as np
import os
import subprocess


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)

class TestSolverWrapperKratosStructure6_0(KratosUnittest.TestCase):

    def test_solver_wrapper_kratos_structure_6_0_tube3d_static(self):

        # Create set up files for kratos
        print_box('setup Kratos static case')

        dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)), f'test_kratos_structure_6_0_tube3d_static')
        set_up_file = os.path.join(dir_tmp, 'set_up_kratos.sh')
        p = subprocess.Popen(f'sh {set_up_file}', cwd=dir_tmp, shell=True)
        p.wait()


        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_kratos_structure_6_0_tube3d_static', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        solver_param = parameters['solver_wrappers'][0]

        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            solver_param['settings'].SetString('working_directory',
                                               f'test_kratos_structure_6_0_tube3d_static/kratos_structure')

        print_box('started tests for Kratos Tube3D static')


        solver = CreateInstance(solver_param)

        load_interface = solver.GetInterfaceInput()

        load_interface_list=  load_interface.GetPythonList()
        self.nr_of_nodes = int(len(load_interface_list) / 4)
        self.load_list = []

        initial_pressure = 10000

        pressure = initial_pressure
        traction_x = 100.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])
        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)

        solver.Initialize()

        solver.InitializeSolutionStep()
        output_1 = solver.SolveSolutionStep(load_interface).deepcopy()
        solver.FinalizeSolutionStep()

        pressure *= 0.5
        traction *= 0.5

        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)
        solver.InitializeSolutionStep()
        output_2 = solver.SolveSolutionStep(load_interface).deepcopy()
        solver.FinalizeSolutionStep()

        pressure = initial_pressure
        traction = np.array([traction_x, traction_y, traction_z])

        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)
        solver.InitializeSolutionStep()
        output_3 = solver.SolveSolutionStep(load_interface).deepcopy()
        solver.FinalizeSolutionStep()

        solver.Finalize()

        # normalize data and compare
        a1 = output_1.GetNumpyArray()
        a2 = output_2.GetNumpyArray()
        a3 = output_3.GetNumpyArray()

        mean = np.mean(a1)
        ref = np.abs(a1 - mean).max()

        a1n = (a1 - mean) / ref
        a2n = (a2 - mean) / ref
        a3n = (a3 - mean) / ref

        self.assertNotAlmostEqual(np.sum(np.abs(a1n - a2n)) / a1n.size, 0., delta=1e-06)
        for i in range(a1.size):
            self.assertAlmostEqual(a1n[i] - a3n[i], 0., delta=1e-07)

    def test_solver_wrapper_kratos_structure_6_0_tube3d_dynamic(self):

        print_box('setup Kratos dynamic case')

        dir_tmp = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                               f'test_kratos_structure_6_0_tube3d_dynamic')
        set_up_file = os.path.join(dir_tmp, 'set_up_kratos.sh')
        p = subprocess.Popen(f'sh {set_up_file}', cwd=dir_tmp, shell=True)
        p.wait()

        # Create set up files for kratos

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_kratos_structure_6_0_tube3d_dynamic', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        solver_param = parameters['solver_wrappers'][0]

        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            solver_param['settings'].SetString('working_directory',
                                               f'test_kratos_structure_6_0_tube3d_dynamic/kratos_structure')

        print_box('started tests for Kratos Tube3D Dynamic')

        solver = CreateInstance(solver_param)

        load_interface = solver.GetInterfaceInput()

        load_interface_list = load_interface.GetPythonList()
        self.nr_of_nodes = int(len(load_interface_list) / 4)

        initial_pressure = 10000

        pressure = initial_pressure
        traction_x = 100.0
        traction_y = 0.0
        traction_z = 0.0
        traction = np.array([traction_x, traction_y, traction_z])
        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)

        solver.Initialize()

        solver.InitializeSolutionStep()
        output_1 = solver.SolveSolutionStep(load_interface).deepcopy()

        pressure *= 0.5
        traction *= 0.5
        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)
        output_2 = solver.SolveSolutionStep(load_interface).deepcopy()

        pressure = initial_pressure
        traction = np.array([traction_x, traction_y, traction_z])
        load_list = self.GetUniformLoadList(pressure, traction)
        load_interface.SetPythonList(load_list)
        output_3 = solver.SolveSolutionStep(load_interface).deepcopy()
        solver.FinalizeSolutionStep()


        solver.Finalize()

        # normalize data and compare
        a1 = output_1.GetNumpyArray()
        a2 = output_2.GetNumpyArray()
        a3 = output_3.GetNumpyArray()

        mean = np.mean(a1)
        ref = np.abs(a1 - mean).max()

        a1n = (a1 - mean) / ref
        a2n = (a2 - mean) / ref
        a3n = (a3 - mean) / ref

        self.assertNotAlmostEqual(np.sum(np.abs(a1n - a2n)) / a1n.size, 0., delta=1e-06)
        for i in range(a1.size):
            self.assertAlmostEqual(a1n[i] - a3n[i], 0., delta=1e-07)

    def GetUniformLoadList(self, pressure, traction):
        load = []
        for i in range(0, self.nr_of_nodes):
            load.append(pressure)
        for i in range(0, self.nr_of_nodes):
            load += traction.tolist()

        return load



if __name__ == '__main__':
    KratosUnittest.main()
