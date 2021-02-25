from coconut import data_structure
from coconut.tools import create_instance

import unittest
import numpy as np
from copy import deepcopy
import os
import math
import multiprocessing
import json
from subprocess import check_call, DEVNULL

import coconut.coupling_components.solver_wrappers.openfoam.open_foam_io as of_io


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


class TestSolverWrapperOpenFoam41(unittest.TestCase):

    def setUp(self):

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_tube_3d', 'test_tube_3d_parameters.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
        self.par_solver = parameters['solver_wrappers'][0]
        self.mp_name_in = self.par_solver['settings']['interface_input'][0]['model_part']
        self.mp_name_out = self.par_solver['settings']['interface_output'][0]['model_part']


        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            self.par_solver['settings']['working_directory'] = 'test_tube_3d'

        self.folder_path = os.path.join(os.getcwd(), self.par_solver['settings']['working_directory'])
        self.delta_t = self.par_solver['settings']['delta_t']
        self.t_prec = self.par_solver['settings']['time_precision']
        self.clean_case()
        self.set_up_case()

    def clean_case(self):
        cmd = os.path.join(self.folder_path, 'Allclean; ')
        cmd += 'rm -rf ' + f'{0:.{self.t_prec}f}; '
        cmd += 'rm -f *.coco'

        check_call(cmd, shell=True, cwd=self.folder_path)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allrun'), shell=True)

    # Test if order of nodes vary with different decomposition
    def test_model_part_nodes_with_different_cores(self):

        print_box("Testing model part nodes with different partitioning")
        # create two solvers with different flow solver partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [1, multiprocessing.cpu_count()]:
            self.par_solver['settings']['cores'] = cores
            solver = create_instance(self.par_solver)
            solver.initialize()
            solver.finalize()
            mp = solver.model.get_model_part(self.mp_name_in)
            x0.append(mp.x0)
            y0.append(mp.y0)
            z0.append(mp.z0)
            ids.append(mp.id)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    # test if nodes are moved to the correct position
    def test_displacement_on_nodes(self):

        print_box("Testing imposed node (radial) displacement")

        # test if nodes are moved to the correct position
        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.par_solver['settings']['cores'] = cores
            solver = create_instance(self.par_solver)
            dr_mag = 0.001
            l0 = 0.05
            mp = solver.model.get_model_part(self.mp_name_in)

            dr = dr_mag * np.sin(2 * np.pi / l0 * mp.x0)
            node_angle = np.arctan2(mp.z0, mp.y0)
            dx = np.zeros(mp.y0.shape)
            dy = dr * np.cos(node_angle)
            dz = dr * np.sin(node_angle)
            displacement = np.column_stack((dx, dy, dz))
            x = mp.x0 + dx
            y = mp.y0 + dy
            z = mp.z0 + dz
            node_coords_ref = np.column_stack((x, y, z))
            solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)

            # update position by iterating once in solver
            solver.initialize()
            solver.initialize_solution_step()
            solver.solve_solution_step(solver.get_interface_input())
            solver.finalize_solution_step()
            solver.finalize()

            if cores > 1:
                check_call(f'cd {self.folder_path} && reconstructPar -latestTime -noFields', shell=True, stdout=DEVNULL)

            node_coords = of_io.get_boundary_points(solver.working_directory, f'{self.delta_t:.{self.t_prec}f}', 'mantle')
            np.testing.assert_allclose(node_coords, node_coords_ref, rtol=1e-12)
            self.clean_case()

    # check if different partitioning gives the same pressure and traction results
    def test_pressure_wall_shear_on_nodes_parallel(self):
        print_box("Testing if pressure/wall shear stress are imposed properly when run in parallel ")
        output_list = []

        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.par_solver['settings']['cores'] = cores
            solver = create_instance(self.par_solver)
            mp = solver.model.get_model_part(self.mp_name_in)
            displacement = np.zeros((mp.size, 3))
            solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)

            # update position by iterating once in solver
            solver.initialize()
            solver.initialize_solution_step()
            output_interface = solver.solve_solution_step(solver.get_interface_input())
            output_list.append(output_interface.get_interface_data())
            solver.finalize_solution_step()
            solver.finalize()
            self.clean_case()

        ref_output = output_list[0]  # single core result
        max_value = np.max(np.abs(ref_output))
        for output in output_list[1:]:
            np.testing.assert_allclose(output / max_value, ref_output / max_value, atol=1e-10, rtol=0)

    #
    #     #
    #     #         # test if same coordinates always gives same pressure & traction
    #     #         if True:
    #     #             # adapt Parameters, create solver
    #     #             par_solver = deepcopy(par_solver_0)
    #     #             par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
    #     #             par_solver['settings'].SetInt('flow_iterations', 500)
    #     #             solver = cs_tools.CreateInstance(par_solver)
    #     #             solver.Initialize()
    #     #             solver.InitializeSolutionStep()
    #     #
    #     #             # change grid to position 1
    #     #             mp = solver.model['beamoutside_nodes']
    #     #             for node in mp.Nodes:
    #     #                 node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #             output1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
    #     #
    #     #             # change grid to position 2
    #     #             for node in mp.Nodes:
    #     #                 node.Y = 0.005 - 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #             output2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
    #     #
    #     #             # change grid back to position 1
    #     #             for node in mp.Nodes:
    #     #                 node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #             output3 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
    #     #
    #     #             solver.FinalizeSolutionStep()
    #     #             solver.Finalize()
    #     #
    #     #             # normalize data and compare
    #     #             a1 = output1.GetNumpyArray()
    #     #             a2 = output2.GetNumpyArray()
    #     #             a3 = output3.GetNumpyArray()
    #     #
    #     #             mean = np.mean(a1)
    #     #             ref = np.abs(a1 - mean).max()
    #     #
    #     #             a1n = (a1 - mean) / ref
    #     #             a2n = (a2 - mean) / ref
    #     #             a3n = (a3 - mean) / ref
    #     #
    #     #             self.assertNotAlmostEqual(np.sum(np.abs(a1n - a2n)) / a1n.size, 0., delta=1e-12)
    #     #             for i in range(a1.size):
    #     #                 self.assertAlmostEqual(a1n[i] - a3n[i], 0., delta=1e-12)
    #     #
    #     #         # test if correct number of displacements is applied
    #     #         if True:
    #     #             # adapt Parameters, create solver
    #     #             par_solver = deepcopy(par_solver_0)
    #     #             par_solver['settings'].SetInt('flow_iterations', 5)
    #     #             solver = cs_tools.CreateInstance(par_solver)
    #     #
    #     #             # give value to DISPLACEMENT variable
    #     #             mp = solver.model['beamoutside_nodes']
    #     #             for node in mp.Nodes:
    #     #                 dy = 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
    #     #
    #     #             # run solver for some timesteps and iterations
    #     #             solver.Initialize()
    #     #             timesteps = 3
    #     #             iterations = 4
    #     #             for i in range(timesteps):
    #     #                 solver.InitializeSolutionStep()
    #     #                 for j in range(iterations):
    #     #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
    #     #                 solver.FinalizeSolutionStep()
    #     #             solver.Finalize()
    #     #
    #     #             # create solver to check coordinates at last timestep
    #     #             par_solver['settings'].SetInt('timestep_start', timesteps)
    #     #             solver = cs_tools.CreateInstance(par_solver)
    #     #             solver.Initialize()
    #     #             solver.Finalize()
    #     #
    #     #             # check if displacement was applied correct number of times
    #     #             mp = solver.model['beamoutside_nodes']
    #     #             for node in mp.Nodes:
    #     #                 n = timesteps * iterations
    #     #                 y_goal = 0.005 + n * 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #                 self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)
    #     #
    #     #         # test if restart option works correctly
    #     #         if True:
    #     #             # adapt Parameters, create solver
    #     #             par_solver = deepcopy(par_solver_0)
    #     #             par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
    #     #             par_solver['settings'].SetInt('flow_iterations', 500)
    #     #             solver = cs_tools.CreateInstance(par_solver)
    #     #
    #     #             # give value to DISPLACEMENT variable
    #     #             mp = solver.model['beamoutside_nodes']
    #     #             for node in mp.Nodes:
    #     #                 dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
    #     #
    #     #             # run solver for 2 timesteps
    #     #             solver.Initialize()
    #     #             for i in range(4):
    #     #                 solver.InitializeSolutionStep()
    #     #                 for j in range(2):
    #     #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
    #     #                 solver.FinalizeSolutionStep()
    #     #             solver.Finalize()
    #     #
    #     #             # get data for solver without restart
    #     #             interface1 = solver.GetInterfaceOutput().deepcopy()
    #     #             data1 = interface1.GetNumpyArray().copy()
    #     #
    #     #             # create solver which restarts at timestep 2
    #     #             par_solver['settings'].SetInt('timestep_start', 2)
    #     #             solver = cs_tools.CreateInstance(par_solver)
    #     #
    #     #             # give value to DISPLACEMENT variable
    #     #             mp = solver.model['beamoutside_nodes']
    #     #             for node in mp.Nodes:
    #     #                 dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
    #     #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
    #     #
    #     #             # run solver for 2 more timesteps
    #     #             solver.Initialize()
    #     #             for i in range(2):
    #     #                 solver.InitializeSolutionStep()
    #     #                 for j in range(2):
    #     #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
    #     #                 solver.FinalizeSolutionStep()
    #     #             solver.Finalize()
    #     #
    #     #             # get data for solver with restart
    #     #             interface2 = solver.GetInterfaceOutput().deepcopy()
    #     #             data2 = interface2.GetNumpyArray().copy()
    #     #
    #     #             # compare coordinates of Nodes
    #     #             mp1 = interface1.model['beamoutside_nodes']
    #     #             mp2 = interface2.model['beamoutside_nodes']
    #     #             for node1 in mp1.Nodes:
    #     #                 node2 = mp2.GetNode(node1.Id)
    #     #                 self.assertAlmostEqual(node1.X, node2.X, delta=1e-16)
    #     #                 self.assertAlmostEqual(node1.Y, node2.Y, delta=1e-16)
    #     #                 self.assertAlmostEqual(node1.Z, node2.Z, delta=1e-16)
    #     #
    #     #             # normalize pressure and traction data and compare
    #     #             mean = np.mean(data1)
    #     #             ref = np.abs(data1 - mean).max()
    #     #
    #     #             data1n = (data1 - mean) / ref
    #     #             data2n = (data2 - mean) / ref
    #     #
    #     #             for i in range(data1n.size):
    #     #                 # print(f'data1n = {data1n[i]:.1e}, data2n = {data2n[i]:.1e}, diff = {data1n[i] - data2n[i]:.1e}')
    #     #                 self.assertAlmostEqual(data1n[i] - data2n[i], 0., delta=1e-14)
    #
    #     print('Finishing tests for OpenFOAM/4.1 tube_3d.')


if __name__ == '__main__':
    unittest.main()  # If this script is executed directly (and NOT through "test_CoSimulationApplication.py') the "KratosUnitTest.py"-script is launched which executes all functions in the class.
