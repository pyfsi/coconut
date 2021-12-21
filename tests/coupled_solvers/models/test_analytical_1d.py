from coconut import tools
from coconut.tools import create_instance

import unittest
import numpy as np
import json
import os
import shutil


class TestAnalytical1D(unittest.TestCase):

    def setUp(self):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to coupled_solvers directory

        # read settings
        parameter_file_name = os.path.join(dir_name, 'test_analytical_1d.json')
        with open(parameter_file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.solver_models_parameters = self.parameters['settings']['solver_models']
        self.model_part_name_input = self.solver_models_parameters[0]['settings']['interface_input'][0]['model_part']
        self.variable_input = self.solver_models_parameters[0]['settings']['interface_input'][0]['variables'][0]
        self.model_part_name_output = self.solver_models_parameters[1]['settings']['interface_output'][0]['model_part']
        self.variable_output = self.solver_models_parameters[1]['settings']['interface_output'][0]['variables'][0]

        # set working directories
        self.working_dir = os.path.join(dir_name, 'coupled_solver_tmp')
        working_dir_cfd = os.path.join(self.working_dir, 'CFD')
        working_dir_csm = os.path.join(self.working_dir, 'CSM')
        self.solver_models_parameters[0]['settings']['working_directory'] = os.path.relpath(working_dir_cfd,
                                                                                            start=self.working_dir)
        self.solver_models_parameters[1]['settings']['working_directory'] = os.path.relpath(working_dir_csm,
                                                                                            start=self.working_dir)

        # setup
        shutil.rmtree(os.path.join(dir_name, self.working_dir), ignore_errors=True)
        os.mkdir(self.working_dir)
        os.mkdir(working_dir_cfd)
        os.mkdir(working_dir_csm)
        shutil.copy(os.path.join(dir_name, '../setup_tube_flow/solver_parameters.json'), working_dir_cfd)
        shutil.copy(os.path.join(dir_name, '../setup_tube_structure/solver_parameters.json'), working_dir_csm)

        # create solvers
        for i in range(2):
            for key in ('delta_t', 'timestep_start'):
                self.solver_models_parameters[i]['settings'][key] = self.parameters['settings'][key]
        with tools.cd(self.working_dir):
            self.flow_solver = tools.create_instance(self.solver_models_parameters[0])  # flow model
            self.structure_solver = tools.create_instance(self.solver_models_parameters[1])  # structure model

        self.interface_input = self.flow_solver.get_interface_input()  # residual of displacements
        self.interface_output = self.structure_solver.get_interface_output()  # displacements
        self.m = self.interface_input.size // 3

    def test_1d_analytical(self):
        for key in ('delta_t', 'timestep_start'):
            tools.remove_recursively(key, self.solver_models_parameters)
        with tools.cd(self.working_dir):
            model = create_instance(self.parameters)

        model.initialize()
        model.initialize_solution_step()

        # test is_ready
        self.assertTrue(model.is_ready())

        # test filter_q
        input = np.zeros((self.m, 3))
        input[:, 1] = np.random.rand(self.m)
        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, input)
        np.testing.assert_array_almost_equal(model.filter_q(self.interface_input).get_interface_data(),
                                             np.zeros(self.m * 3))

        # test jacobian
        for solver in (self.flow_solver, self.structure_solver):
            solver.initialize()
            solver.initialize_solution_step()

        # zero deformation
        x = model.x.get_variable_data(self.model_part_name_input, self.variable_input)
        dx = np.zeros_like(x)
        dx[:, 1] = 1e-10

        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, x - dx / 2)
        intermediate = self.flow_solver.solve_solution_step(self.interface_input)
        xt1 = self.structure_solver.solve_solution_step(intermediate).get_interface_data()

        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, x + dx / 2)
        intermediate = self.flow_solver.solve_solution_step(self.interface_input)
        xt2 = self.structure_solver.solve_solution_step(intermediate).get_interface_data()

        dxt = xt2 - xt1
        dr = dxt - dx.flatten()

        self.interface_input.set_interface_data(dr)
        prediction_dxt = model.predict(self.interface_input).get_interface_data()

        self.assertLess(np.linalg.norm(prediction_dxt - dxt), 1e-8)

        # random deformation
        x = np.zeros((self.m, 3))
        x[:, 1] = np.random.rand(self.m) * self.flow_solver.d * 0.01
        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, x)
        self.interface_output.set_variable_data(self.model_part_name_output, self.variable_output, x)

        model.add(self.interface_input * 0, self.interface_output)
        model.finalize_solution_step()
        model.initialize_solution_step()

        self.flow_solver.solve_solution_step(self.interface_input)
        intermediate = self.flow_solver.solve_solution_step(self.interface_input)
        self.structure_solver.solve_solution_step(intermediate)
        for solver in (self.flow_solver, self.structure_solver):
            solver.finalize_solution_step()
            solver.initialize_solution_step()

        dx = np.zeros_like(x)
        dx[:, 1] = 1e-10

        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, x - dx / 2)
        intermediate = self.flow_solver.solve_solution_step(self.interface_input)
        xt1 = self.structure_solver.solve_solution_step(intermediate).get_interface_data()

        self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, x + dx / 2)
        intermediate = self.flow_solver.solve_solution_step(self.interface_input)
        xt2 = self.structure_solver.solve_solution_step(intermediate).get_interface_data()

        dxt = xt2 - xt1
        dr = dxt - dx.flatten()

        self.interface_input.set_interface_data(dr)
        prediction_dxt = model.predict(self.interface_input).get_interface_data()

        self.assertLess(np.linalg.norm(prediction_dxt - dxt), 1e-8)

    def test_update_iteration(self):
        for key in ('delta_t', 'timestep_start'):
            tools.remove_recursively(key, self.solver_models_parameters)
        with tools.cd(self.working_dir):
            model_no_update = create_instance(self.parameters)
        self.parameters['settings']['update_every_iteration'] = True
        for key in ('delta_t', 'timestep_start'):
            tools.remove_recursively(key, self.solver_models_parameters)
        with tools.cd(self.working_dir):
            model_update = create_instance(self.parameters)

        for model, check in zip((model_no_update, model_update), (True, False)):
            model.initialize()
            model.initialize_solution_step()

            input = np.zeros((self.m, 3))
            input[:, 1] = np.random.rand(self.m)

            self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, input)
            prediction_1 = model.predict(self.interface_input)

            r_in = np.zeros((self.m, 3))
            r_in[:, 1] = np.random.rand(self.m)
            self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, r_in)
            xt_in = np.zeros((self.m, 3))
            xt_in[:, 1] = np.random.rand(self.m)
            self.interface_output.set_variable_data(self.model_part_name_output, self.variable_output, xt_in)
            model.add(self.interface_input, self.interface_output)

            self.interface_input.set_variable_data(self.model_part_name_input, self.variable_input, input)
            prediction_2 = model.predict(self.interface_input)

            self.assertEqual(prediction_1 == prediction_2, check)

    def tearDown(self):
        shutil.rmtree(self.working_dir)

if __name__ == '__main__':
    unittest.main()
