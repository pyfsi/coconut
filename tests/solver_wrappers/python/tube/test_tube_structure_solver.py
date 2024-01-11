from coconut.tools import create_instance
import coconut.coupling_components.solver_wrappers.python.banded as bnd

import unittest
import os
import shutil
import json
import numpy as np


def get_dp(x):
    return 1500 * np.sin(2 * np.pi / 0.05 * x)


class TestSolverWrapperTubeStructureSolver(unittest.TestCase):

    def setUp(self):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to python.tube directory

        # read settings
        parameter_file_name = os.path.join(dir_name, 'test_tube_structure/test_tube_structure_solver.json')
        with open(parameter_file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.variable = 'pressure'
        self.model_part_name = 'wall'

        # set working directory
        self.working_dir = os.path.join(dir_name, 'test_tube_structure/CSM')
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)

        # setup
        shutil.rmtree(os.path.join(dir_name, self.working_dir), ignore_errors=True)
        os.mkdir(self.working_dir)
        shutil.copy(os.path.join(dir_name, 'test_tube_structure/setup_tube_structure/solver_parameters.json'),
                    self.working_dir)

    def test_pressure(self):
        # test if same pressure always give same displacement

        # create solver
        solver_1 = create_instance(self.parameters)
        solver_2 = create_instance(self.parameters)
        solvers = [solver_1, solver_2]
        for solver in solvers:
            solver.initialize()
            solver.initialize_solution_step()

        # change solver_1 to end pressure and solve
        interface_input = solver_1.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = get_dp(model_part.z0).reshape(-1, 1)
        interface_input.set_variable_data(self.model_part_name, self.variable, data)
        output1_end = solver_1.solve_solution_step(interface_input).copy()

        # change solver_2 to intermediate pressure and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = 0.5 * get_dp(model_part.z0).reshape(-1, 1)
        interface_input.set_variable_data(self.model_part_name, self.variable, data)
        solver_2.solve_solution_step(interface_input).copy()

        # change solver_2 to end position and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = get_dp(model_part.z0).reshape(-1, 1)
        interface_input.set_variable_data(self.model_part_name, self.variable, data)
        output2_end = solver_2.solve_solution_step(interface_input).copy()

        for solver in solvers:
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

        # compare
        a1 = output1_end.get_interface_data()
        a2 = output2_end.get_interface_data()

        np.testing.assert_allclose(a1, a2, atol=1e-15)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['save_restart'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set pressure
        z0 = interface_input.get_model_part(self.model_part_name).z0
        pressure = interface_input.get_variable_data(self.model_part_name, 'pressure')
        pressure[:, 0] = get_dp(z0)

        # run solver for 4 time steps
        for i in range(4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.model_part_name, 'pressure', i * pressure)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        displacement_1 = interface_output.get_variable_data(self.model_part_name, 'displacement')

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # run solver for 2 more time steps
        for i in range(2, 4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.model_part_name, 'pressure', i * pressure)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        solver.finalize()

        # get data for solver with restart
        interface_output = solver.get_interface_output()
        displacement_2 = interface_output.get_variable_data(self.model_part_name, 'displacement')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if pressure and traction are equal
        np.testing.assert_allclose(displacement_1, displacement_2, rtol=1e-14)

    def test_coupling_convergence(self):
        # test if check of coupling convergence works correctly

        # adapt parameters, create solver
        solver = create_instance(self.parameters)
        solver.check_coupling_convergence = True
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set pressure
        z0 = interface_input.get_model_part(self.model_part_name).z0
        pressure = interface_input.get_variable_data(self.model_part_name, 'pressure')
        pressure[:, 0] = get_dp(z0)

        solver.initialize_solution_step()

        # coupling iteration
        interface_input.set_variable_data(self.model_part_name, 'pressure', pressure)
        solver.solve_solution_step(interface_input)

        self.assertTrue(solver.coupling_convergence)

        solver.output_solution_step()
        solver.finalize_solution_step()
        solver.finalize()

    def test_jacobian(self):
        # create solver and set pressure
        solver = self.set_up_solver()

        # verify derivative
        r = np.array(solver.r)
        dr = 1e-3 * np.ones_like(r)
        dr[:2] = 0
        dr[-2:] = 0
        solver.r = r - dr / 2
        f1 = solver.get_residual()
        solver.r = r + dr / 2
        f2 = solver.get_residual()
        solver.r = r
        j = bnd.to_dense(solver.get_jacobian())
        ju = j[2:-2, 2:-2]
        f1 = f1[2:-2]
        f2 = f2[2:-2]
        dr = dr[2:-2]
        d1 = f2 - f1
        d2 = ju @ dr
        self.assertLess(np.linalg.norm(d1 - d2), 1e-16)

    def test_surrogate_jacobian(self):
        # test surrogate Jacobian

        # create solver and set displacement
        solver = self.set_up_solver()

        # verify surrogate Jacobian
        input = np.zeros((solver.m, 1))
        p = np.array(solver.p)
        dp = 1e5 * np.ones_like(p)

        input[:, 0] = p - dp / 2
        interface_input = solver.get_interface_input()
        interface_input.set_variable_data(self.model_part_name, self.variable, input)
        f1 = solver.solve_solution_step(interface_input).get_variable_data(self.model_part_name, 'displacement')[:, 1]

        input[:, 0] = p + dp / 2
        interface_input = solver.get_interface_input()
        interface_input.set_variable_data(self.model_part_name, self.variable, input)
        f2 = solver.solve_solution_step(interface_input).get_variable_data(self.model_part_name, 'displacement')[:, 1]

        jsurrogate = solver.get_surrogate_jacobian()
        d1 = f2 - f1
        d2 = jsurrogate @ dp
        print(np.linalg.norm(d1 - d2))
        self.assertLess(np.linalg.norm(d1 - d2), 1e-16)

    def set_up_solver(self):
        # create solver
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()

        # set pressure
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[get_dp(model_part.z0[i])]for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        solver.solve_solution_step(interface_input)

        return solver

    def tearDown(self):
        shutil.rmtree(self.working_dir)


if __name__ == '__main__':
    unittest.main()
