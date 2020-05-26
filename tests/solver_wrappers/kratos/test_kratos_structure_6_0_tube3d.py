from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance
from coconut.coupling_components.interface import Interface
from coconut.coupling_components import output_vtk_process

import numpy as np
import os


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)

class TestSolverWrapperKratosStructure6_0(KratosUnittest.TestCase):

    def test_solver_wrapper_kratos_structure_6_0_tube3d_static(self):

        print_box('started tests for Kratos Tube3D static')

        # Create set up files for kratos
        os.system("cd test_kratos_structure_6_0_tube3d_static; sh set_up_kratos.sh")

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_kratos_structure_6_0_tube3d_static', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        solver_param = parameters['solver_wrappers'][0]
        solver = CreateInstance(solver_param)

        load_interface = solver.GetInterfaceInput()

        load_interface_list=  load_interface.GetPythonList()
        self.nr_of_nodes = int(len(load_interface_list) / 4)
        self.load_list = []

        initial_pressure = 10000

        self.pressure = initial_pressure
        traction_x = 100.0
        traction_y = 0.0
        traction_z = 0.0
        self.traction = np.array([traction_x, traction_y, traction_z])
        self.ApplyUniformLoad()
        load_interface.SetPythonList(self.load_list)

        solver.Initialize()

        out_obj = output_vtk_process.OutputInterfaceProcess(load_interface,'test_kratos_structure_6_0_tube3d_static/kratos_structure')

        solver.InitializeSolutionStep()
        disp_out_1 = solver.SolveSolutionStep(load_interface).GetNumpyArray()
        solver.FinalizeSolutionStep()
        out_obj.Write()

        self.pressure *= 0.5
        self.traction *= 0.5

        self.ApplyUniformLoad()
        load_interface.SetPythonList(self.load_list)
        solver.InitializeSolutionStep()
        solver.SolveSolutionStep(load_interface)
        solver.FinalizeSolutionStep()
        out_obj.Write()

        self.pressure = initial_pressure
        self.traction = np.array([traction_x, traction_y, traction_z])

        self.ApplyUniformLoad()
        load_interface.SetPythonList(self.load_list)
        solver.InitializeSolutionStep()
        disp_out_2 = solver.SolveSolutionStep(load_interface).GetNumpyArray()
        solver.FinalizeSolutionStep()
        out_obj.Write()

        solver.Finalize()
        os.system("mv *.msh test_kratos_structure_6_0_tube3d_static/kratos_structure; mv *.res test_kratos_structure_6_0_tube3d_static/kratos_structure; mv *.lst  test_kratos_structure_6_0_tube3d_static/kratos_structure")
        #os.system("cd test_kratos_structure_6_0_tube3d_static/kratos_structure; python3 convert_gid_to_vtk.py; cd -")

        norm_disp_out_1 = np.linalg.norm(disp_out_1.size)
        #norm_disp_diff = np.linalg.norm(disp_out_1 - disp_out_2)
        # print('Relative norm difference: ', norm_disp_diff / norm_disp_out_1)
        # print('Absolute norm difference: ', norm_disp_diff)

        for i in range(disp_out_1.size):
            self.assertAlmostEqual((disp_out_1[i] - disp_out_2[i])/norm_disp_out_1, 0., delta=1e-12)

    def test_solver_wrapper_kratos_structure_6_0_tube3d_dynamic(self):

        print_box('started tests for Kratos Tube3D Dynamic')

        # Create set up files for kratos
        os.system("cd test_kratos_structure_6_0_tube3d_dynamic; sh set_up_kratos.sh")

        parameter_file_name = os.path.join(os.path.dirname(__file__),
                                           'test_kratos_structure_6_0_tube3d_dynamic', 'test_solver_wrapper.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        solver_param = parameters['solver_wrappers'][0]
        solver = CreateInstance(solver_param)

        load_interface = solver.GetInterfaceInput()

        load_interface_list = load_interface.GetPythonList()
        self.nr_of_nodes = int(len(load_interface_list) / 4)
        self.load_list = []

        initial_pressure = 10000

        self.pressure = initial_pressure
        traction_x = 100.0
        traction_y = 0.0
        traction_z = 0.0
        self.traction = np.array([traction_x, traction_y, traction_z])
        self.ApplyUniformLoad()
        load_interface.SetPythonList(self.load_list)

        solver.Initialize()

        solver.InitializeSolutionStep()
        disp_out_1 = solver.SolveSolutionStep(load_interface).GetNumpyArray()
        #solver.FinalizeSolutionStep()


        self.ApplyUniformLoad()
        load_interface.SetPythonList(self.load_list)
        #solver.InitializeSolutionStep()
        disp_out_2 = solver.SolveSolutionStep(load_interface).GetNumpyArray()
        solver.FinalizeSolutionStep()

        solver.Finalize()
        os.system("mv *.msh test_kratos_structure_6_0_tube3d_dynamic/kratos_structure; mv *.res test_kratos_structure_6_0_tube3d_dynamic/kratos_structure; mv *.lst  test_kratos_structure_6_0_tube3d_dynamic/kratos_structure")
        os.system("cd test_kratos_structure_6_0_tube3d_dynamic/kratos_structure; python3 convert_gid_to_vtk.py; cd -")

        norm_disp_out_1 = np.linalg.norm(disp_out_1.size)
        # norm_disp_diff = np.linalg.norm(disp_out_1 - disp_out_2)
        # print('Relative norm difference: ', norm_disp_diff / norm_disp_out_1)
        # print('Absolute norm difference: ', norm_disp_diff)

        for i in range(disp_out_1.size):
            self.assertAlmostEqual((disp_out_1[i] - disp_out_2[i]) / norm_disp_out_1, 0., delta=1e-12)

    def ApplyUniformLoad(self):
        self.load_list = []
        for i in range(0, self.nr_of_nodes):
            self.load_list.append(self.pressure)
        for i in range(0, self.nr_of_nodes):
            self.load_list +=  self.traction.tolist()



if __name__ == '__main__':
    test_obj = TestSolverWrapperKratosStructure6_0()
    test_obj.test_solver_wrapper_kratos_structure_6_0_tube3d_static()
    test_obj.test_solver_wrapper_kratos_structure_6_0_tube3d_dynamic()


