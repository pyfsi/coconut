from coconut.tools import create_instance
import coconut.coupling_components.solver_wrappers.python.banded as bnd

import unittest
import os
import shutil
import json
import numpy as np


def get_dy(x):
    return 0.0001 * np.sin(2 * np.pi / 0.05 * x)


class TestSolverWrapperTubeFlowSolver(unittest.TestCase):

    def setUp(self):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to python.tube directory

        # read settings
        parameter_file_name = os.path.join(dir_name, 'test_tube_flow/test_tube_flow_solver.json')
        with open(parameter_file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.variable = 'displacement'
        self.model_part_name = 'wall'

        # set working directory
        self.working_dir = os.path.join(dir_name, 'test_tube_flow/CFD')
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)

        # setup
        shutil.rmtree(os.path.join(dir_name, self.working_dir), ignore_errors=True)
        os.mkdir(self.working_dir)
        shutil.copy(os.path.join(dir_name, 'test_tube_flow/setup_tube_flow/solver_parameters.json'), self.working_dir)

    def test_pressure(self):
        # test if same coordinates always give same pressure

        # create solver
        solver_1 = create_instance(self.parameters)
        solver_2 = create_instance(self.parameters)
        solvers = [solver_1, solver_2]
        for solver in solvers:
            solver.initialize()
            solver.initialize_solution_step()

        # change solver_1 to end position and solve
        interface_input = solver_1.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[0, get_dy(model_part.z0[i]), 0] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        output1_end = solver_1.solve_solution_step(interface_input).copy()

        # change solver_2 to intermediate position and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[0, -get_dy(model_part.z0[i]), 0] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        solver_2.solve_solution_step(interface_input).copy()

        # change solver_2 to end position and solve
        interface_input = solver_2.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[0, get_dy(model_part.z0[i]), 0] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        output2_end = solver_2.solve_solution_step(interface_input).copy()

        for solver in solvers:
            solver.finalize_solution_step()
            solver.output_solution_step()
            solver.finalize()

        # compare
        a1 = output1_end.get_interface_data()
        a2 = output2_end.get_interface_data()

        np.testing.assert_allclose(a1, a2, atol=1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['save_restart'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        z0 = interface_input.get_model_part(self.model_part_name).z0
        displacement = interface_input.get_variable_data(self.model_part_name, 'displacement')
        displacement[:, 1] = get_dy(z0)

        # run solver for 4 time steps
        for i in range(4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.model_part_name, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        pressure_1 = interface_output.get_variable_data(self.model_part_name, 'pressure')
        traction_1 = interface_output.get_variable_data(self.model_part_name, 'traction')

        # create solver which restarts at time step 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        solver.initialize()
        interface_input = solver.get_interface_input()

        # run solver for 2 more time steps
        for i in range(2, 4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.model_part_name, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        solver.finalize()

        # get data for solver with restart
        interface_output = solver.get_interface_output()
        pressure_2 = interface_output.get_variable_data(self.model_part_name, 'pressure')
        traction_2 = interface_output.get_variable_data(self.model_part_name, 'traction')

        # check if undeformed coordinate (coordinates of model part) are equal
        self.assertTrue(interface_x_1.has_same_model_parts(interface_x_2))
        self.assertTrue(interface_y_1.has_same_model_parts(interface_y_2))

        # check if pressure and traction are equal
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-14)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-14)

    def test_coupling_convergence(self):
        # test if check of coupling convergence works correctly

        # adapt parameters, create solver
        self.parameters['settings']['newtonmax'] = 5
        solver = create_instance(self.parameters)
        solver.check_coupling_convergence = True
        solver.initialize()
        interface_input = solver.get_interface_input()

        # set displacement
        z0 = interface_input.get_model_part(self.model_part_name).z0
        displacement = interface_input.get_variable_data(self.model_part_name, 'displacement')
        displacement[:, 1] = get_dy(z0)

        solver.initialize_solution_step()

        # first coupling iteration
        interface_input.set_variable_data(self.model_part_name, 'displacement', displacement)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # second coupling iteration
        interface_input.set_variable_data(self.model_part_name, 'displacement', 2 * displacement)
        solver.solve_solution_step(interface_input)

        self.assertFalse(solver.coupling_convergence)

        # third coupling iteration
        interface_input.set_variable_data(self.model_part_name, 'displacement', 2 * displacement)
        solver.solve_solution_step(interface_input)

        self.assertTrue(solver.coupling_convergence)

        solver.output_solution_step()
        solver.finalize_solution_step()
        solver.finalize()

    def test_area_jacobian(self):
        # test area derivative

        # create solver and set displacement
        solver = self.set_up_solver()

        # verify derivative
        a = np.array(solver.a)
        da = 1e-8 * np.ones_like(a)
        da[0] = 0
        da[-1] = 0
        solver.a = a - da / 2
        f1 = solver.get_residual()
        solver.a = a + da / 2
        f2 = solver.get_residual()
        solver.a = a
        ja = solver.get_jacobian_area()
        f1 = f1[2:-2]
        f2 = f2[2:-2]
        da = da[1:-1]
        d1 = f2 - f1
        d2 = ja @ da
        self.assertLess(np.linalg.norm(d1 - d2), 1e-15)

    def test_pressure_jacobian(self):
        # test pressure derivative

        # create solver and set displacement
        solver = self.set_up_solver()

        # verify derivative
        p = np.array(solver.p)
        dp = 1e-5 * np.ones_like(p)
        dp[0] = 0
        dp[-1] = 0
        solver.p = p - dp / 2
        f1 = solver.get_residual()
        solver.p = p + dp / 2
        f2 = solver.get_residual()
        solver.p = p
        j = bnd.to_dense(solver.get_jacobian())
        jp = j[2:-2, 3:-2:2]
        f1 = f1[2:-2]
        f2 = f2[2:-2]
        dp = dp[1:-1]
        d1 = f2 - f1
        d2 = jp @ dp
        self.assertLess(np.linalg.norm(d1 - d2), 1e-15)

    def test_velocity_jacobian(self):
        # test velocity derivative

        # create solver and set displacement
        solver = self.set_up_solver()

        # verify derivative
        u = np.array(solver.u)
        du = 1e-5 * np.ones_like(u)
        du[0] = 0
        du[-1] = 0
        solver.u = u - du / 2
        f1 = solver.get_residual()
        solver.u = u + du / 2
        f2 = solver.get_residual()
        solver.u = u
        j = bnd.to_dense(solver.get_jacobian())
        ju = j[2:-2, 2:-2:2]
        f1 = f1[2:-2]
        f2 = f2[2:-2]
        du = du[1:-1]
        d1 = f2 - f1
        d2 = ju @ du
        self.assertLess(np.linalg.norm(d1 - d2), 1e-15)

    def test_surrogate_jacobian(self):
        # test surrogate Jacobian

        # create solver and set displacement
        solver = self.set_up_solver()

        # verify surrogate Jacobian
        input_data = np.zeros((solver.m, 3))
        a = np.array(solver.a[1:-1])
        r = np.sqrt(a / np.pi)  # local radius
        dr = 1e-18 * np.ones_like(r)

        input_data[:, 1] = r - dr / 2 - solver.d / 2
        interface_input = solver.get_interface_input()
        interface_input.set_variable_data(self.model_part_name, self.variable, input_data)
        f1 = solver.solve_solution_step(interface_input).get_variable_data(self.model_part_name, 'pressure').flatten()

        input_data[:, 1] = r + dr / 2 - solver.d / 2
        interface_input = solver.get_interface_input()
        interface_input.set_variable_data(self.model_part_name, self.variable, input_data)
        f2 = solver.solve_solution_step(interface_input).get_variable_data(self.model_part_name, 'pressure').flatten()

        jsurrogate = solver.get_surrogate_jacobian()
        d1 = f2 - f1
        d2 = jsurrogate @ dr
        self.assertLess(np.linalg.norm(d1 - d2), 5e-7)

    def set_up_solver(self):
        # create solver
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()

        # set displacement
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.model_part_name)
        data = [[0, get_dy(model_part.z0[i]), 0] for i in range(model_part.size)]
        interface_input.set_variable_data(self.model_part_name, self.variable, np.array(data))
        solver.solve_solution_step(interface_input)

        return solver

    def tearDown(self):
        shutil.rmtree(self.working_dir)


if __name__ == '__main__':
    unittest.main()
