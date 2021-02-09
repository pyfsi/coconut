from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance

import numpy as np
from copy import deepcopy
import os
import math
import multiprocessing
from subprocess import check_call, DEVNULL

import coconut.coupling_components.solver_wrappers.openfoam.open_foam_io as of_io


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


class TestSolverWrapperOpenFoam41(KratosUnittest.TestCase):

    def setUp(self):

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_tube_3d', 'test_tube_3d_parameters.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        self.par_solver = parameters['solver_wrappers'][0]


        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            self.par_solver['settings'].SetString('working_directory', 'test_tube_3d')

        self.displacement = vars(data_structure)['DISPLACEMENT']

        self.folder_path = os.path.join(self.par_solver['settings']['working_directory'].GetString())
        self.dt = self.par_solver['settings']['dt'].GetDouble()
        self.t_prec = self.par_solver['settings']['time_precision'].GetDouble()
        self.clean_case()
        self.set_up_case()

    def clean_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allclean'), shell=True)
        check_call('rm -rf ' + os.path.join(self.folder_path, f'{0:.{self.t_prec}f}'), shell=True)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allrun'), shell=True)

    # Test if order of nodes vary with different decomposition
    def test_model_part_nodes_with_different_cores(self):

        print_box("Testing model part nodes with different partitioning")
        # create two solvers with different flow solver partitioning
        model_parts = []
        for cores in [1, multiprocessing.cpu_count()]:
            self.par_solver['settings'].SetInt('cores', cores)
            solver = CreateInstance(self.par_solver)
            solver.Initialize()
            solver.Finalize()
            model_parts.append(deepcopy(solver.model['mantle_input']))

        # compare Nodes in ModelParts between both solvers
        mp1, mp2 = model_parts
        for node1 in mp1.Nodes:
            node2 = mp2.GetNode(node1.Id)
            self.assertEqual(node1.X0, node2.X0)
            self.assertEqual(node1.Y0, node2.Y0)
            self.assertEqual(node1.Z0, node2.Z0)

    # test if nodes are moved to the correct position
    def test_displacement_on_nodes(self):

        print_box("Testing imposed node (radial) displacement")

        # test if nodes are moved to the correct position
        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.par_solver['settings'].SetInt('cores', cores)
            solver = CreateInstance(self.par_solver)
            mp = solver.model['mantle_input']
            dr_mag = 0.001
            r0 = 0.005
            l0 = 0.05

            for node in mp.Nodes:
                dr = dr_mag * np.sin(2 * np.pi / l0 * node.X0)
                node_angle = np.arctan2(node.Z0, node.Y0)
                dy = dr * np.cos(node_angle)
                dz = dr * np.sin(node_angle)
                node.SetSolutionStepValue(self.displacement, 0, [0., dy, dz])


            # update position by iterating once in solver
            solver.Initialize()
            solver.InitializeSolutionStep()
            solver.SolveSolutionStep(solver.GetInterfaceInput())
            solver.FinalizeSolutionStep()
            solver.Finalize()

            if cores > 1:
                check_call(f'cd {self.folder_path} && reconstructPar -latestTime -noFields', shell=True, stdout=DEVNULL)


            node_coords = of_io.get_boundary_points(solver.working_directory, f'{self.dt:.{self.t_prec}f}', 'mantle')

            for i, node in enumerate(mp.Nodes):
                r_goal = r0 + dr_mag * np.sin(2 * np.pi / l0 * node.X0)
                y = node_coords[i, 1]
                z = node_coords[i, 2]
                self.assertAlmostEqual(math.sqrt(y ** 2 + z ** 2), r_goal, delta=1e-15)

            self.clean_case()

    # # check if different partitioning gives the same pressure and traction results
    def test_pressure_wall_shear_on_nodes_parallel(self):
        print_box("Testing if perssure/wall shear stress are imposed properly when run in parallel ")
        output_list = []

        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.par_solver['settings'].SetInt('cores', cores)
            solver = CreateInstance(self.par_solver)

            mp = solver.model['mantle_input']

            for node in mp.Nodes:
                node.SetSolutionStepValue(self.displacement, 0, [0., 0., 0.])

            # update position by iterating once in solver
            solver.Initialize()
            solver.InitializeSolutionStep()
            output = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
            output_list.append(output.GetNumpyArray())
            solver.FinalizeSolutionStep()
            solver.Finalize()
            self.clean_case()

        ref_output = output_list[0]  # single core result
        for output in output_list[1:]:
            np.testing.assert_array_almost_equal(output / np.max(np.abs(ref_output)),
                                                 ref_output / np.max(np.abs(ref_output)), decimal=6)

        #
        #         # test if same coordinates always gives same pressure & traction
        #         if True:
        #             # adapt Parameters, create solver
        #             par_solver = deepcopy(par_solver_0)
        #             par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
        #             par_solver['settings'].SetInt('flow_iterations', 500)
        #             solver = cs_tools.CreateInstance(par_solver)
        #             solver.Initialize()
        #             solver.InitializeSolutionStep()
        #
        #             # change grid to position 1
        #             mp = solver.model['beamoutside_nodes']
        #             for node in mp.Nodes:
        #                 node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
        #             output1 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
        #
        #             # change grid to position 2
        #             for node in mp.Nodes:
        #                 node.Y = 0.005 - 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
        #             output2 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
        #
        #             # change grid back to position 1
        #             for node in mp.Nodes:
        #                 node.Y = 0.005 + 0.0005 * np.sin(2 * np.pi / 0.05 * node.X)
        #             output3 = solver.SolveSolutionStep(solver.GetInterfaceInput()).deepcopy()
        #
        #             solver.FinalizeSolutionStep()
        #             solver.Finalize()
        #
        #             # normalize data and compare
        #             a1 = output1.GetNumpyArray()
        #             a2 = output2.GetNumpyArray()
        #             a3 = output3.GetNumpyArray()
        #
        #             mean = np.mean(a1)
        #             ref = np.abs(a1 - mean).max()
        #
        #             a1n = (a1 - mean) / ref
        #             a2n = (a2 - mean) / ref
        #             a3n = (a3 - mean) / ref
        #
        #             self.assertNotAlmostEqual(np.sum(np.abs(a1n - a2n)) / a1n.size, 0., delta=1e-12)
        #             for i in range(a1.size):
        #                 self.assertAlmostEqual(a1n[i] - a3n[i], 0., delta=1e-12)
        #
        #         # test if correct number of displacements is applied
        #         if True:
        #             # adapt Parameters, create solver
        #             par_solver = deepcopy(par_solver_0)
        #             par_solver['settings'].SetInt('flow_iterations', 5)
        #             solver = cs_tools.CreateInstance(par_solver)
        #
        #             # give value to DISPLACEMENT variable
        #             mp = solver.model['beamoutside_nodes']
        #             for node in mp.Nodes:
        #                 dy = 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
        #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
        #
        #             # run solver for some timesteps and iterations
        #             solver.Initialize()
        #             timesteps = 3
        #             iterations = 4
        #             for i in range(timesteps):
        #                 solver.InitializeSolutionStep()
        #                 for j in range(iterations):
        #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
        #                 solver.FinalizeSolutionStep()
        #             solver.Finalize()
        #
        #             # create solver to check coordinates at last timestep
        #             par_solver['settings'].SetInt('timestep_start', timesteps)
        #             solver = cs_tools.CreateInstance(par_solver)
        #             solver.Initialize()
        #             solver.Finalize()
        #
        #             # check if displacement was applied correct number of times
        #             mp = solver.model['beamoutside_nodes']
        #             for node in mp.Nodes:
        #                 n = timesteps * iterations
        #                 y_goal = 0.005 + n * 0.0001 * np.sin(2 * np.pi / 0.05 * node.X)
        #                 self.assertAlmostEqual(node.Y, y_goal, delta=1e-16)
        #
        #         # test if restart option works correctly
        #         if True:
        #             # adapt Parameters, create solver
        #             par_solver = deepcopy(par_solver_0)
        #             par_solver['settings'].SetInt('cores', multiprocessing.cpu_count())
        #             par_solver['settings'].SetInt('flow_iterations', 500)
        #             solver = cs_tools.CreateInstance(par_solver)
        #
        #             # give value to DISPLACEMENT variable
        #             mp = solver.model['beamoutside_nodes']
        #             for node in mp.Nodes:
        #                 dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
        #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
        #
        #             # run solver for 2 timesteps
        #             solver.Initialize()
        #             for i in range(4):
        #                 solver.InitializeSolutionStep()
        #                 for j in range(2):
        #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
        #                 solver.FinalizeSolutionStep()
        #             solver.Finalize()
        #
        #             # get data for solver without restart
        #             interface1 = solver.GetInterfaceOutput().deepcopy()
        #             data1 = interface1.GetNumpyArray().copy()
        #
        #             # create solver which restarts at timestep 2
        #             par_solver['settings'].SetInt('timestep_start', 2)
        #             solver = cs_tools.CreateInstance(par_solver)
        #
        #             # give value to DISPLACEMENT variable
        #             mp = solver.model['beamoutside_nodes']
        #             for node in mp.Nodes:
        #                 dy = 0.0002 * np.sin(2 * np.pi / 0.05 * node.X)
        #                 node.SetSolutionStepValue(displacement, 0, [0., dy, 0.])
        #
        #             # run solver for 2 more timesteps
        #             solver.Initialize()
        #             for i in range(2):
        #                 solver.InitializeSolutionStep()
        #                 for j in range(2):
        #                     solver.SolveSolutionStep(solver.GetInterfaceInput())
        #                 solver.FinalizeSolutionStep()
        #             solver.Finalize()
        #
        #             # get data for solver with restart
        #             interface2 = solver.GetInterfaceOutput().deepcopy()
        #             data2 = interface2.GetNumpyArray().copy()
        #
        #             # compare coordinates of Nodes
        #             mp1 = interface1.model['beamoutside_nodes']
        #             mp2 = interface2.model['beamoutside_nodes']
        #             for node1 in mp1.Nodes:
        #                 node2 = mp2.GetNode(node1.Id)
        #                 self.assertAlmostEqual(node1.X, node2.X, delta=1e-16)
        #                 self.assertAlmostEqual(node1.Y, node2.Y, delta=1e-16)
        #                 self.assertAlmostEqual(node1.Z, node2.Z, delta=1e-16)
        #
        #             # normalize pressure and traction data and compare
        #             mean = np.mean(data1)
        #             ref = np.abs(data1 - mean).max()
        #
        #             data1n = (data1 - mean) / ref
        #             data2n = (data2 - mean) / ref
        #
        #             for i in range(data1n.size):
        #                 # print(f'data1n = {data1n[i]:.1e}, data2n = {data2n[i]:.1e}, diff = {data1n[i] - data2n[i]:.1e}')
        #                 self.assertAlmostEqual(data1n[i] - data2n[i], 0., delta=1e-14)

        print('Finishing tests for OpenFOAM/4.1 tube_3d.')


if __name__ == '__main__':
    KratosUnittest.main()  # If this script is executed directly (and NOT through "test_CoSimulationApplication.py') the "KratosUnitTest.py"-script is launched which executes all functions in the class.
