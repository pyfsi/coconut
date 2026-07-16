import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json
import multiprocessing
import shutil

from coconut.tools import create_instance, get_solver_env, rm_timed


# ==============================================================================
# 1. FLUENT SOLID WRAPPER (Stefan 2-Phase Case)
# ==============================================================================
class TestSolverWrapperPCSolidFluent(unittest.TestCase):
    version = 'xxxxRx'
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, f'test_v{cls.version}/stefan_2phase/solid/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/stefan_2phase/solid/CFD')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/stefan_2phase/solid/setup_fluent')

        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)
            fluent_solver_module = f'fluent.v{cls.version}'
            env = get_solver_env(fluent_solver_module, cls.working_dir)
            # FIX: Execute explicitly with bash to avoid permission 126 errors
            subprocess.check_call('bash setup_fluent.sh', shell=True, cwd=cls.working_dir, env=env)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = self.working_dir
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in_faces = 'boundary_in_faces'
        self.mp_name_in_nodes = 'boundary_in_nodes'
        self.mp_name_out_nodes = 'boundary_out_nodes'

    @classmethod
    def tearDownClass(cls):
        pass
        # if cls.setup_case:
        #     rm_timed(cls.working_dir)

    def test_fluent_melting_displacement(self):
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        interface_input = solver.get_interface_input()
        heat_flux = interface_input.get_variable_data(self.mp_name_in_faces, 'heat_flux')
        heat_flux[:] = -50000.0
        interface_input.set_variable_data(self.mp_name_in_faces, 'heat_flux', heat_flux)

        solver.initialize_solution_step()
        interface_output = solver.solve_solution_step(interface_input)

        disp = interface_output.get_variable_data(self.mp_name_out_nodes, 'displacement')
        max_disp = np.max(np.linalg.norm(disp, axis=1))
        self.assertTrue(max_disp > 0.0)

        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()


# ==============================================================================
# 2. FLUENT LIQUID WRAPPER (No RB) (Stefan 2-Phase Case)
# ==============================================================================
class TestSolverWrapperPCLiquid(unittest.TestCase):
    version = 'xxxxRx'
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, f'test_v{cls.version}/stefan_2phase/liquid/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/stefan_2phase/liquid/CFD')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/stefan_2phase/liquid/setup_fluent')

        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)
            fluent_solver_module = f'fluent.v{cls.version}'
            env = get_solver_env(fluent_solver_module, cls.working_dir)
            # FIX: Execute explicitly with bash to avoid permission 126 errors
            subprocess.check_call('bash setup_fluent.sh', shell=True, cwd=cls.working_dir, env=env)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = self.working_dir
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in = 'boundary_in_nodes'
        self.mp_name_out_faces = 'boundary_out_faces'

    @classmethod
    def tearDownClass(cls):
        pass
        # if cls.setup_case:
        #     rm_timed(cls.working_dir)

    def test_liquid_interface_tracking(self):
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        interface_input = solver.get_interface_input()
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 0] += 1e-5
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        solver.initialize_solution_step()
        interface_output = solver.solve_solution_step(interface_input)

        heat_flux = interface_output.get_variable_data(self.mp_name_out_faces, 'heat_flux')
        self.assertEqual(heat_flux.shape[1], 1)

        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()


# ==============================================================================
# 3. FLUENT LIQUID WRAPPER WITH RIGID BODY (Hirata Case)
# ==============================================================================
class TestSolverWrapperPCLiquidRB(unittest.TestCase):
    version = 'xxxxRx'
    setup_case = True

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))
        cls.file_name = join(dir_name, f'test_v{cls.version}/hirata/liquid/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/hirata/liquid/CFD')
        cls.setup_dir = join(dir_name, f'test_v{cls.version}/hirata/liquid/setup_fluent')

        if cls.setup_case:
            shutil.rmtree(cls.working_dir, ignore_errors=True)
            shutil.copytree(cls.setup_dir, cls.working_dir)
            fluent_solver_module = f'fluent.v{cls.version}'
            env = get_solver_env(fluent_solver_module, cls.working_dir)
            # FIX: Execute explicitly with bash to avoid permission 126 errors
            subprocess.check_call('bash setup_fluent.sh', shell=True, cwd=cls.working_dir, env=env)

    def setUp(self):
        with open(self.file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = self.working_dir
        self.parameters['settings']['cores'] = min(4, multiprocessing.cpu_count())
        self.mp_name_in = 'inner_in_nodes'
        self.mp_name_out_faces = 'inner_out_faces'

    @classmethod
    def tearDownClass(cls):
        pass
        # if cls.setup_case:
        #     rm_timed(cls.working_dir)

    def test_data_transfer_and_rbm(self):
        self.parameters['settings']['flow_iterations'] = 1
        self.parameters['settings']['RB']['iteration_max'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        interface_input = solver.get_interface_input()
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] -= 1e-6
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        solver.initialize_solution_step()
        interface_output = solver.solve_solution_step(interface_input)
        heat_flux = interface_output.get_variable_data(self.mp_name_out_faces, 'heat_flux')
        self.assertEqual(heat_flux.shape[1], 1)

        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

    def test_contact_model(self):
        self.parameters['settings']['flow_iterations'] = 1
        self.parameters['settings']['RB']['iteration_max'] = 1
        solver = create_instance(self.parameters)
        solver.initialize()

        interface_input = solver.get_interface_input()
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] -= 0.05
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        solver.initialize_solution_step()
        solver.solve_solution_step(interface_input)

        self.assertTrue(len(solver.rbm_solver.contact_patches) > 0)
        self.assertTrue(solver.rbm_solver.h_min < solver.rbm_solver.h_ul)

        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

    def test_fictitious_mass(self):
        self.parameters['settings']['flow_iterations'] = 1
        self.parameters['settings']['RB']['iteration_max'] = 1
        self.parameters['settings']['RB']['fict_coeff'] = 8.0

        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()

        interface_input = solver.get_interface_input()
        solver.solve_solution_step(interface_input)

        m_sys_active = solver.rbm_solver.M_sys
        self.assertTrue(m_sys_active > 0.0)

        solver.finalize_solution_step()
        solver.output_solution_step()
        solver.finalize()

    def test_restart(self):
        # --- DEEL 1: CONTINUOUS RUN ---
        self.parameters['settings']['flow_iterations'] = 2
        self.parameters['settings']['RB']['iteration_max'] = 2

        # 1. Forceer de solver om na elke tijdstap een restart-file op te slaan
        self.parameters['settings']['save_restart'] = 1

        solver = create_instance(self.parameters)
        solver.initialize()

        interface_input = solver.get_interface_input()
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        for i in range(3):
            solver.initialize_solution_step()
            displacement[:, 1] -= 1e-6
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
            solver.output_solution_step()

        hf_1 = solver.get_interface_output().get_variable_data(self.mp_name_out_faces, 'heat_flux')
        coords_1 = solver.get_coordinates()[self.mp_name_in]['coords']
        solver.finalize()

        # --- DEEL 2: RESTART RUN ---
        self.parameters['settings']['timestep_start'] = 2

        # 2. Vertel Fluent expliciet welk bestand hij moet inlezen
        self.parameters['settings']['case_file'] = 'case_timestep2.cas.h5'

        # 3. Vertel de Rigid Body physics solver dat hij de pickle array moet laden
        self.parameters['settings']['RB']['restart_rb'] = 2

        solver_restart = create_instance(self.parameters)
        solver_restart.initialize()
        interface_input_res = solver_restart.get_interface_input()

        for i in range(2, 3):
            solver_restart.initialize_solution_step()
            # De verplaatsing vector moet nog een extra stapje ondergaan ten opzichte van start
            displacement[:, 1] -= 1e-6
            interface_input_res.set_variable_data(self.mp_name_in, 'displacement', displacement)
            solver_restart.solve_solution_step(interface_input_res)
            solver_restart.finalize_solution_step()
            solver_restart.output_solution_step()

        hf_2 = solver_restart.get_interface_output().get_variable_data(self.mp_name_out_faces, 'heat_flux')
        coords_2 = solver_restart.get_coordinates()[self.mp_name_in]['coords']
        solver_restart.finalize()

        # --- DEEL 3: VERGELIJKING ---
        np.testing.assert_allclose(coords_1, coords_2, rtol=1e-5, atol=1e-9, err_msg="Restart coordinates mismatch!")
        np.testing.assert_allclose(hf_1, hf_2, rtol=1e-5, atol=1e-9, err_msg="Restart heat flux mismatch!")