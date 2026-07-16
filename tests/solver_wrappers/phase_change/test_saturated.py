import unittest
import numpy as np
import os
from os.path import join
import json
import shutil

from coconut.tools import create_instance, rm_timed


class TestSolverWrapperPCSolidPython(unittest.TestCase):
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, 'test_v2024R2/hirata/solid/parameters.json')
        cls.working_dir = join(dir_name, 'test_v2024R2/hirata/solid/CFD')
        cls.setup_dir = join(dir_name, 'test_v2024R2/hirata/solid')

        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            os.makedirs(cls.working_dir, exist_ok=True)

            for file in ['solver_parameters.json', 'interface_nodes.csv']:
                src = join(cls.setup_dir, file)
                if os.path.exists(src):
                    shutil.copy(src, cls.working_dir)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)

        # FIX: Gebruik het absolute pad zodat de Python solver de mappen altijd vindt
        self.parameters['settings']['working_directory'] = self.working_dir
        self.mp_name_in_faces = 'inner_in_faces'
        self.mp_name_in_nodes = 'inner_in_nodes'
        self.mp_name_out_nodes = 'inner_out_nodes'

    @classmethod
    def tearDownClass(cls):
        pass
        # if cls.setup_case:
        #     rm_timed(cls.working_dir)

    def test_melting_displacement(self):
        solver = create_instance(self.parameters)
        solver.initialize()

        q_val = -10000.0

        interface_input = solver.get_interface_input()
        heat_flux = interface_input.get_variable_data(self.mp_name_in_faces, 'heat_flux')
        heat_flux[:] = q_val
        interface_input.set_variable_data(self.mp_name_in_faces, 'heat_flux', heat_flux)

        solver.initialize_solution_step()
        interface_output = solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()
        solver.output_solution_step()

        disp = interface_output.get_variable_data(self.mp_name_out_nodes, 'displacement')

        max_disp = np.max(np.linalg.norm(disp, axis=1))
        self.assertTrue(max_disp > 0.0)

        solver.finalize()

    def test_restart(self):
        """Test if history states (prev_disp) are flawlessly recovered during a restart."""
        # Ensure the mapper has enough iterations to converge fully,
        # preventing microscopic differences between cold and warm starts.
        if 'conservative_mapper' in self.parameters['settings']:
            self.parameters['settings']['conservative_mapper']['mapper_iterations'] = 100

        self.parameters['settings']['save_restart'] = 1
        self.parameters['settings']['timestep_start'] = 0

        solver = create_instance(self.parameters)
        solver.save_restart = 1
        solver.initialize()

        interface_input = solver.get_interface_input()
        heat_flux = interface_input.get_variable_data(self.mp_name_in_faces, 'heat_flux')

        for i in range(3):
            solver.initialize_solution_step()

            # Feedback loop: in a real coupled simulation, the output displacement becomes the new input.
            if i > 0:
                prev_out_disp = solver.get_interface_output().get_variable_data(self.mp_name_out_nodes, 'displacement')
                interface_input.set_variable_data(self.mp_name_in_nodes, 'displacement', prev_out_disp)

            heat_flux[:] = -5000.0
            interface_input.set_variable_data(self.mp_name_in_faces, 'heat_flux', heat_flux)

            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()

        disp_1 = solver.get_interface_output().get_variable_data(self.mp_name_out_nodes, 'displacement')
        solver.finalize()

        pickle_path = join(self.working_dir, 'case_timestep2.pickle')
        self.assertTrue(os.path.exists(pickle_path), f"Restart file was not written to location: {pickle_path}")

        self.parameters['settings']['timestep_start'] = 2

        solver_restart = create_instance(self.parameters)
        solver_restart.save_restart = 1
        solver_restart.initialize()

        interface_input_res = solver_restart.get_interface_input()
        for i in range(2, 3):
            solver_restart.initialize_solution_step()

            # For the first step after restart (i=2), feedback is not needed,
            # as the solver's __init__ already populates interface_input correctly.
            if i > 2:
                prev_out_disp = solver_restart.get_interface_output().get_variable_data(self.mp_name_out_nodes,
                                                                                        'displacement')
                interface_input_res.set_variable_data(self.mp_name_in_nodes, 'displacement', prev_out_disp)

            heat_flux[:] = -5000.0
            interface_input_res.set_variable_data(self.mp_name_in_faces, 'heat_flux', heat_flux)

            solver_restart.solve_solution_step(interface_input_res)
            solver_restart.finalize_solution_step()
            solver_restart.output_solution_step()

        disp_2 = solver_restart.get_interface_output().get_variable_data(self.mp_name_out_nodes, 'displacement')
        solver_restart.finalize()

        diff = np.max(np.linalg.norm(disp_1 - disp_2, axis=1))
        np.testing.assert_allclose(disp_1, disp_2, rtol=1e-5, atol=1e-9,
                                   err_msg=f"Restart displacements do not match! Max diff: {diff}")