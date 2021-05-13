from coconut.tools import create_instance, get_solver_env, solver_available, print_box

import unittest
import numpy as np
import os
import multiprocessing
import json
from subprocess import check_call, DEVNULL

import coconut.coupling_components.solver_wrappers.openfoam.open_foam_io as of_io


version = '41'


@unittest.skipUnless(solver_available(f'openfoam.v{version}'), f'openfoam.v{version} not available')
class TestSolverWrapperOpenFoam41(unittest.TestCase):

    def setUp(self):

        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_tube_3d', 'parameters.json')

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

        solver_name = self.par_solver['type'].replace('solver_wrappers.', '')
        self.env = get_solver_env(solver_name, self.folder_path)
        self.clean_case()
        self.set_up_case()

    def clean_case(self):
        cmd = 'sh ' + os.path.join(self.folder_path, 'Allclean; ')
        cmd += 'rm -rf ' + f'{0:.{self.t_prec}f}; '
        cmd += 'rm -f *.coco'

        check_call(cmd, shell=True, cwd=self.folder_path, env=self.env)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allrun'), shell=True, env=self.env)

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
                check_call(f'cd {self.folder_path} && reconstructPar -latestTime -noFields', shell=True, stdout=DEVNULL, env=self.env)

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



if __name__ == '__main__':
    unittest.main()  # If this script is executed directly (and NOT through "test_CoSimulationApplication.py') the "KratosUnitTest.py"-script is launched which executes all functions in the class.
