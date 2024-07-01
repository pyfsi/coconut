from coconut.tools import create_instance, get_solver_env
import coconut.coupling_components.solver_wrappers.openfoam.openfoam_io as of_io

import unittest
import numpy as np
import os
from os.path import join
import multiprocessing
import shutil
import re
import json
import glob
from subprocess import check_call, DEVNULL


class TestSolverWrapperOpenFOAM(unittest.TestCase):
    version = None  # OpenFOAM version without dot, e.g. 41 , set in sub-class

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to openfoam directory
        cls.file_name = join(dir_name, f'test_v{cls.version}/tube3d/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/tube3d/CFD')

        # setup
        shutil.rmtree(cls.working_dir, ignore_errors=True)
        shutil.copytree(join(dir_name, f'test_v{cls.version}/tube3d/setup_openfoam'), cls.working_dir)

    def setUp(self):
        with open(self.file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory

        self.mp_name_in = self.parameters['settings']['interface_input'][0]['model_part']
        self.mp_name_out = self.parameters['settings']['interface_output'][0]['model_part']

        self.folder_path = os.path.join(os.getcwd(), self.parameters['settings']['working_directory'])
        self.delta_t = self.parameters['settings']['delta_t']
        self.t_prec = self.parameters['settings']['time_precision']
        self.max_cores = min(4, multiprocessing.cpu_count())  # number of cores for parallel calculation

        solver_name = self.parameters['type'].replace('solver_wrappers.', '')
        self.env = get_solver_env(solver_name, self.folder_path)
        self.clean_case()
        self.set_up_case()

    def tearDown(self):
        self.clean_case()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.working_dir)

    def clean_case(self):  # means redoing blockMesh which takes a long time
        check_call('sh ' + os.path.join(self.folder_path, 'Allclean'), shell=True, env=self.env)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'prepareCase'), shell=True, env=self.env)

    def set_cores(self, cores):
        if cores == 1:
            self.parameters['settings']['parallel'] = False
        else:
            self.parameters['settings']['parallel'] = True
            decompose_file_name = os.path.join(self.folder_path, 'system', 'decomposeParDict')
            with open(decompose_file_name, 'r') as f:
                old_dict = f.read()
            new_dict = re.sub(r'numberOfSubdomains[\s\n]+' + of_io.int_pattern, f'numberOfSubdomains   {cores}',
                              old_dict)
            with open(decompose_file_name, 'w') as f:
                f.write(new_dict)

    def set_tolerance(self, p_atol, U_atol):
        solution_file_name = os.path.join(self.folder_path, 'system', 'fvSolution')
        new_tol = f'outerCorrectorResidualControl\n' \
                  f'    {{\n' \
                  f'    p\n' \
                  f'        {{\n' \
                  f'            tolerance   {p_atol:g};\n' \
                  f'            relTol      0;\n' \
                  f'        }}\n' \
                  f'    U\n' \
                  f'        {{\n' \
                  f'            tolerance   {U_atol:g};\n' \
                  f'            relTol      0;\n' \
                  f'        }}\n' \
                  f'    }}'
        with open(solution_file_name, 'r') as f:
            old_dict = f.read()
        new_dict = re.sub(of_io.get_nested_dict(old_dict, 'outerCorrectorResidualControl'), new_tol, old_dict)
        with open(solution_file_name, 'w') as f:
            f.write(new_dict)

    # noinspection PyMethodMayBeStatic
    def get_dy_dz(self, x, y, z):
        dr = 0.0001 * np.sin(2 * np.pi / 0.05 * x)
        theta = np.arctan2(z, y)
        dy = dr * np.cos(theta)
        dz = dr * np.sin(theta)
        return dy, dz

    # noinspection PyMethodMayBeStatic
    def rm_time_dirs(self, solver):
        # remove 0.0000 type folder created by functionObject during serial run,
        # resulting in an error when decomposed for a parallel run in the same directory
        for time_dir in glob.glob(join(solver.working_directory, '0.*')):
            shutil.rmtree(time_dir)

    # test if order of nodes vary with different decomposition
    def test_model_part_nodes_with_different_cores(self):
        # create two solvers with different flow solver partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [1, self.max_cores]:
            self.set_cores(cores)
            solver = create_instance(self.parameters)
            solver.initialize()
            solver.finalize()
            model_part = solver.model.get_model_part(self.mp_name_in)
            x0.append(model_part.x0)
            y0.append(model_part.y0)
            z0.append(model_part.z0)
            ids.append(model_part.id)

            self.rm_time_dirs(solver)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    # test if nodes are moved to the correct position
    def test_displacement_on_nodes(self):
        for cores in [1, self.max_cores]:
            self.set_up_case()
            self.set_cores(cores)
            solver = create_instance(self.parameters)
            solver.initialize()

            model_part = solver.model.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dx = np.zeros(x0.shape)
            dy, dz = self.get_dy_dz(x0, y0, z0)
            displacement = np.column_stack((dx, dy, dz))
            x = x0 + dx
            y = y0 + dy
            z = z0 + dz
            node_coords_ref = np.column_stack((x, y, z))
            interface_disp = solver.get_interface_input()
            interface_disp.set_variable_data(self.mp_name_in, 'displacement', displacement)

            # update position by iterating once in solver
            solver.initialize_solution_step()
            solver.solve_solution_step(interface_disp)
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

            if cores > 1:
                check_call(f'cd {self.folder_path} && reconstructPar -latestTime -noFields', shell=True, stdout=DEVNULL,
                           env=self.env)

            _, node_coords = of_io.get_boundary_points(solver.working_directory,
                                                       f'{self.delta_t:.{solver.time_precision}f}', 'mantle')
            np.testing.assert_allclose(node_coords, node_coords_ref, rtol=1e-12)

            self.rm_time_dirs(solver)

    # check if different partitioning gives the same pressure and traction results
    def test_pressure_wall_shear_on_nodes_parallel(self):
        output_list = []
        for cores in [1, self.max_cores]:
            self.set_up_case()
            self.set_cores(cores)
            solver = create_instance(self.parameters)
            solver.initialize()
            interface_input = solver.get_interface_input()

            model_part = interface_input.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dx = np.zeros(x0.shape)
            dy, dz = self.get_dy_dz(x0, y0, z0)
            displacement = np.column_stack((dx, dy, dz))
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

            # update position by iterating once in solver
            solver.initialize_solution_step()
            output_interface = solver.solve_solution_step(interface_input)
            output_list.append(output_interface.get_interface_data())
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

            self.rm_time_dirs(solver)

        ref_output = output_list[0]  # single core result
        max_value = np.max(np.abs(ref_output))
        for output in output_list[1:]:
            np.testing.assert_allclose(output / max_value, ref_output / max_value, atol=1e-10, rtol=0)

    # test if same coordinates always gives same pressure & traction
    def test_pressure_wall_shear(self):
        # adapt parameters, create solver
        self.set_cores(1)
        solver = create_instance(self.parameters)
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
        solver.output_solution_step()
        solver.finalize()

        pr_amp = 0.5 * (np.max((pressure[0] + pressure[2]) / 2) - np.min((pressure[0] + pressure[2]) / 2))
        tr_amp = 0.5 * (np.max((traction[0] + traction[2]) / 2) - np.min((traction[0] + traction[2]) / 2))

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0] / pr_amp, pressure[2] / pr_amp, atol=1e-4, rtol=0)
        np.testing.assert_allclose(traction[0] / tr_amp, traction[2] / tr_amp, atol=1e-4, rtol=0)

        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-4)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-4)

    # test if restart option works correctly
    def test_restart(self):
        # adapt parameters, create solver
        cores = 4
        self.set_cores(cores)
        self.parameters['settings']['cores'] = cores
        self.parameters['settings']['save_restart'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
        dx = np.zeros(x0.shape)
        displacement = np.column_stack((dx, dy, dz))
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
        nr_time_steps = 4

        # run solver for 4 time steps
        for i in range(nr_time_steps):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()

        if cores > 1:
            check_call(f'reconstructPar -latestTime -noFields', shell=True, cwd=self.folder_path, stdout=DEVNULL,
                       env=self.env)
        _, coords_1 = of_io.get_boundary_points(solver.working_directory,
                                                f'{nr_time_steps * self.delta_t:.{solver.time_precision}f}', 'mantle')
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_1 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # run solver for 2 more time steps
        for i in range(2, nr_time_steps):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        if cores > 1:
            check_call(f'reconstructPar -latestTime -noFields', shell=True, cwd=self.folder_path, stdout=DEVNULL,
                       env=self.env)
        _, coords_2 = of_io.get_boundary_points(solver.working_directory,
                                                f'{nr_time_steps * self.delta_t:.{solver.time_precision}f}', 'mantle')
        solver.finalize()

        # get data for solver with restart
        interface_output = solver.get_interface_output()
        pressure_2 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_2 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if coordinates of ModelParts are equal
        # ==>  check if deformed coordinates are equal
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-15)

        # check if pressure and traction are equal
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-9)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-9)

    def test_coupling_convergence(self):
        # test if check of coupling convergence works correctly

        # change of tolerance is necessary since the test does not succeed with strict tolerances
        # this is because the mesh deformation calculation is slightly nonlinear due to the explicit non-orthogonality
        # correction in the Laplacian calculation, see documentation for more information
        self.set_tolerance(1e-3, 1e-3)

        for cores in (1, 4):
            # adapt parameters, create solver
            self.set_cores(cores)
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
            solver.check_coupling_convergence = True
            solver.initialize()
            interface_input = solver.get_interface_input()

            # set displacement
            model_part = interface_input.get_model_part(self.mp_name_in)
            x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
            dy, dz = self.get_dy_dz(x0, y0, z0)
            dx = np.zeros(x0.shape)
            displacement = np.column_stack((dx, dy, dz))

            solver.initialize_solution_step()

            # first coupling iteration
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            solver.solve_solution_step(interface_input)
            self.assertFalse(solver.coupling_convergence)

            # second coupling iteration
            interface_input.set_variable_data(self.mp_name_in, 'displacement', 2 * displacement)
            solver.solve_solution_step(interface_input)
            self.assertFalse(solver.coupling_convergence)

            # third coupling iteration
            interface_input.set_variable_data(self.mp_name_in, 'displacement', 2 * displacement)
            solver.solve_solution_step(interface_input)
            self.assertTrue(solver.coupling_convergence)

            solver.output_solution_step()
            solver.finalize_solution_step()
            solver.finalize()

            self.rm_time_dirs(solver)


if __name__ == '__main__':
    unittest.main()
