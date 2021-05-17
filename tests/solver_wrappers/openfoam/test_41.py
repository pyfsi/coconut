from coconut.tools import create_instance
from coconut.tools import get_solver_env

import unittest
import os
import multiprocessing
import re
import json
from subprocess import check_call, DEVNULL
import numpy as np

import coconut.coupling_components.solver_wrappers.openfoam.openfoam_io as of_io


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)


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

    def tearDown(self):
        check_call(f"""ps aux | awk '{{if (($0 !~ /awk/) && ($0 ~ /coconut_pimpleFoam/)) system("kill " $2)}}'""", shell=True)

    def clean_case(self):

        check_call('sh ' + os.path.join(self.folder_path, 'Allclean'), shell=True, env=self.env)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allrun'), shell=True, env=self.env)

    def set_cores(self, cores):
        if cores == 1:
            self.par_solver['settings']['parallel'] = False
        else:
            self.par_solver['settings']['parallel'] = True
            decompose_file_name = os.path.join(self.folder_path, 'system', 'decomposeParDict')
            with open(decompose_file_name, 'r') as f:
                dict = f.read()
            new_dict = re.sub(r'numberOfSubdomains[\s\n]+' + of_io.int_pattern, f'numberOfSubdomains   {cores}', dict)

            with open(decompose_file_name, 'w') as f:
                f.write(new_dict)

    def get_dy_dz(self, x, y, z):
        dr = 0.0001 * np.sin(2 * np.pi / 0.05 * x)
        theta = np.arctan2(z, y)
        dy = dr * np.cos(theta)
        dz = dr * np.sin(theta)
        return dy, dz

    ## Test if order of nodes vary with different decomposition
    def test_model_part_nodes_with_different_cores(self):

        print_box("Testing model part nodes with different partitioning")
        # create two solvers with different flow solver partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [1, multiprocessing.cpu_count()]:
            self.set_cores(cores)
            solver = create_instance(self.par_solver)
            solver.initialize()
            solver.finalize()
            model_part = solver.model.get_model_part(self.mp_name_in)
            x0.append(model_part.x0)
            y0.append(model_part.y0)
            z0.append(model_part.z0)
            ids.append(model_part.id)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    ## test if nodes are moved to the correct position
    def test_displacement_on_nodes(self):

        print_box("Testing imposed node (radial) displacement")

        # test if nodes are moved to the correct position
        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.set_cores(cores)
            solver = create_instance(self.par_solver)

            model_part = solver.model.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dx = np.zeros(x0.shape)
            dy, dz = self.get_dy_dz(x0, y0, z0)
            displacement = np.column_stack((dx, dy, dz))
            x = x0 + dx
            y = y0 + dy
            z = z0 + dz
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

            _, node_coords = of_io.get_boundary_points(solver.working_directory,
                                                       f'{self.delta_t:.{solver.time_precision}f}', 'mantle')
            np.testing.assert_allclose(node_coords, node_coords_ref, rtol=1e-12)
            self.clean_case()

    ## check if different partitioning gives the same pressure and traction results
    def test_pressure_wall_shear_on_nodes_parallel(self):
        print_box("Testing if pressure/wall shear stress are imposed properly when run in parallel ")
        output_list = []

        for cores in [1, multiprocessing.cpu_count()]:
            self.set_up_case()
            self.set_cores(cores)
            solver = create_instance(self.par_solver)
            model_part = solver.model.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dx = np.zeros(x0.shape)
            dy, dz = self.get_dy_dz(x0, y0, z0)
            displacement = np.column_stack((dx, dy, dz))
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

    ## test if same coordinates always gives same pressure & traction
    def test_pressure_wall_shear(self):
        print_box("Testing if same displacement gives same pressure/traction")
        # adapt parameters, create solver
        self.set_cores(1)
        solver = create_instance(self.par_solver)
        solver.initialize()
        solver.initialize_solution_step()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
        dydz = np.column_stack((dy, dz))
        dydz_list = [dydz, np.zeros_like(dydz), dydz]

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for dydz in dydz_list:
            displacement[:, 1:] = dydz
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.finalize()

        pr_amp = 0.5*(np.max((pressure[0] + pressure[2])/2) - np.min((pressure[0] + pressure[2])/2))
        tr_amp = 0.5*(np.max((traction[0] + traction[2])/2) - np.min((traction[0] + traction[2])/2))

        # # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0]/pr_amp, pressure[2]/pr_amp, atol=1e-4, rtol=0)
        np.testing.assert_allclose(traction[0]/tr_amp, traction[2]/tr_amp, atol=1e-4, rtol=0)

        #check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-4)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-4)

        #print(np.c_[pressure[0], pressure[1], pressure[2]])

        ## test if restart option works correctly
    def test_restart(self):
        # adapt parameters, create solver
        print_box("Testing restart")
        cores = 4
        self.set_cores(cores)
        self.par_solver['settings']['cores'] = cores
        solver = create_instance(self.par_solver)
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
        dx = np.zeros(x0.shape)
        displacement = np.column_stack((dx, dy, dz))
        solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)
        nr_time_steps = 4
        # run solver for 4 timesteps
        solver.initialize()
        for i in range(nr_time_steps):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()

        if cores > 1:
            check_call(f'reconstructPar -latestTime -noFields', shell=True,cwd=self.folder_path, stdout=DEVNULL, env=self.env)
        _, coords_1 = of_io.get_boundary_points(solver.working_directory,
                                             f'{nr_time_steps * self.delta_t:.{solver.time_precision}f}',
                                             'mantle')
        solver.finalize()
        # get data for solver without restart
        interface_output = solver.get_interface_output()
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_1 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # create solver which restarts at timestep 2
        self.par_solver['settings']['timestep_start'] = 2
        solver = create_instance(self.par_solver)
        interface_input = solver.get_interface_input()

        # run solver for 2 more timesteps
        solver.initialize()
        for i in range(2, nr_time_steps):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        if cores > 1:
            check_call(f'reconstructPar -latestTime -noFields', shell=True,cwd=self.folder_path, stdout=DEVNULL, env=self.env)
        _, coords_2 = of_io.get_boundary_points(solver.working_directory,
                                             f'{nr_time_steps * self.delta_t:.{solver.time_precision}f}',
                                             'mantle')
        solver.finalize()

        # # get data for solver with restart
        interface_output = solver.get_interface_output()
        pressure_2 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_2 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # # check if undeformed coordinate (coordinates of model part) are equal
        # self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        # self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if coordinates of ModelParts are equal
        # ==>  check if deformed coordinates are equal
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-15)

        # # check if pressure and traction are equal
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-9)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-9)




if __name__ == '__main__':
    unittest.main()  # If this script is executed directly (and NOT through "test_CoSimulationApplication.py') the "KratosUnitTest.py"-script is launched which executes all functions in the class.
