from coconut.tools import create_instance

import unittest
import numpy as np
import os
from os.path import join
import subprocess
import json

# TODO: issue: id is a Python build-in function... use ids instead?

# TODO: put somewhere in docs?
"""
- The unittests in this file (i.e. for v2019R1) can be run
directly, using

    python -m unittest -vb solver_wrappers/fluent/test_v2019R1.py

- The unittests for the other Fluent versions currently 
inherit from these classes. Because those modules import
these super classes, they are also run and crash, which is
kinda ugly. 
To run other Fluent version unittests, adapt the file
__ini__.py and use

    python -m unittest discover -vb
    
- The unittests for different Fluent versions can never be
run simultaneously, because the correct software must be 
loaded in advance. 
"""


class TestSolverWrapperFluent2019R1Tube2D(unittest.TestCase):
    version = '2019R1'
    setup_case = True

    @classmethod
    def setUpClass(cls):
        if cls.setup_case:
            dir_tmp = join(os.getcwd(), f'solver_wrappers/fluent/test_v{cls.version}/tube2d')
            p = subprocess.Popen(join(dir_tmp, 'setup_fluent.sh'), cwd=dir_tmp, shell=True)
            p.wait()

    def setUp(self):
        file_name = f'solver_wrappers/fluent/test_v{self.version}/tube2d/parameters.json'
        with open(file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.mp_name_in = 'beamoutside_nodes'
        self.mp_name_out = 'beamoutside_faces'

    def get_dy(self, x):
        return 0.0001 * np.sin(2 * np.pi / 0.05 * x)

    def test_move_nodes(self):
        # test if nodes are moved to the correct position

        # adapt parameters, create solver
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)

        # set displacement
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        x0, y0 = model_part.x0, model_part.y0
        dy = self.get_dy(x0)

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        solver.initialize()
        solver.initialize_solution_step()
        solver.solve_solution_step(interface_input)
        solver.finalize_solution_step()
        coord_data = solver.get_coordinates()
        solver.finalize()

        # check if correct displacement was given
        y = coord_data[self.mp_name_in]['coords'][:, 1]
        np.testing.assert_allclose(y, y0 + dy, rtol=1e-15)

    def test_partitioning(self):
        # test if different partitioning gives the same ModelParts

        # create two solvers with different partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [0, 1]:
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
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

    def test_pressure_traction(self):
        # test if same coordinates always give same pressure & traction

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 0
        self.parameters['settings']['flow_iterations'] = 500
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()
        interface_input = solver.get_interface_input()

        # set displacement
        x0 = interface_input.get_model_part(self.mp_name_in).x0
        dys = [self.get_dy(x0), np.zeros_like(x0), self.get_dy(x0)]
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for dy in dys:
            displacement[:, 1] = dy
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.finalize()

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0], pressure[2], rtol=1e-10)
        np.testing.assert_allclose(traction[0], traction[2], rtol=1e-11)

        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-12)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 30
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # set displacement
        x0 = interface_input.get_model_part(self.mp_name_in).x0
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = self.get_dy(x0)

        # run solver for 4 timesteps
        solver.initialize()
        for i in range(4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        coords_1 = solver.get_coordinates()[self.mp_name_in]['coords']
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_1 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # create solver which restarts at timestep 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # run solver for 2 more timesteps
        solver.initialize()
        for i in range(2, 4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        coords_2 = solver.get_coordinates()[self.mp_name_in]['coords']
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
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-14)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-14)


class TestSolverWrapperFluent2019R1Tube3D(unittest.TestCase):
    version = '2019R1'
    setup_case = True

    @classmethod
    def setUpClass(cls):
        if cls.setup_case:
            dir_tmp = join(os.getcwd(),f'solver_wrappers/fluent/test_v{cls.version}/tube3d')
            p = subprocess.Popen(join(dir_tmp, 'setup_fluent.sh'), cwd=dir_tmp, shell=True)
            p.wait()

    def setUp(self):
        file_name = f'solver_wrappers/fluent/test_v{self.version}/tube3d/parameters.json'
        with open(file_name) as parameter_file:
            self.parameters = json.load(parameter_file)
        self.mp_name_in = 'wall_nodes'
        self.mp_name_out = 'wall_faces'

    def get_dy_dz(self, x, y, z):
        dr = 0.0001 * np.sin(2 * np.pi / 0.05 * x)
        theta = np.arctan2(z, y)
        dy = dr * np.cos(theta)
        dz = dr * np.sin(theta)
        return dy, dz

    def test_move_nodes(self):
        # test if nodes are moved to the correct position

        # adapt parameters, create solver
        self.parameters['settings']['flow_iterations'] = 1
        solver = create_instance(self.parameters)

        # set displacement
        interface_input = solver.get_interface_input()
        model_part = interface_input.get_model_part(self.mp_name_in)

        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        displacement[:, 2] = dz
        interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        solver.initialize()
        solver.initialize_solution_step()
        solver.solve_solution_step(solver.get_interface_input())
        solver.finalize_solution_step()
        coord_data = solver.get_coordinates()
        solver.finalize()

        # check if correct displacement was given
        y = coord_data[self.mp_name_in]['coords'][:, 1]
        z = coord_data[self.mp_name_in]['coords'][:, 2]
        np.testing.assert_allclose(y, y0 + dy, rtol=1e-15)
        np.testing.assert_allclose(z, z0 + dz, rtol=1e-15)

    def test_partitioning(self):
        # test if different partitioning gives the same ModelParts

        # create two solvers with different partitioning
        x0, y0, z0, ids = [], [], [], []
        for cores in [0, 1]:
            self.parameters['settings']['cores'] = cores
            solver = create_instance(self.parameters)
            solver.initialize()
            solver.finalize()
            mp = solver.model.get_model_part('wall_nodes')
            x0.append(mp.x0)
            y0.append(mp.y0)
            z0.append(mp.z0)
            ids.append(mp.id)

        # compare ModelParts of both solvers
        for attr in [x0, y0, z0, ids]:
            np.testing.assert_array_equal(attr[0], attr[1])

    def test_pressure_traction(self):
        # test if same coordinates always give same pressure & traction

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 0
        self.parameters['settings']['flow_iterations'] = 500
        solver = create_instance(self.parameters)
        solver.initialize()
        solver.initialize_solution_step()
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
        dydz = np.column_stack((dy, dz))
        dydzs = [dydz, np.zeros_like(dydz), dydz]
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for dydz in dydzs:
            displacement[:, 1:] = dydz
            interface_input.set_variable_data(self.mp_name_in, 'displacement', displacement)
            interface_output = solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
            traction.append(interface_output.get_variable_data(self.mp_name_out, 'traction'))
        solver.finalize_solution_step()
        solver.finalize()

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0], pressure[2], rtol=1e-11)
        np.testing.assert_allclose(traction[0], traction[2], rtol=1e-11)

        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-12)

        t01 = np.linalg.norm(traction[0] - traction[1])
        t02 = np.linalg.norm(traction[0] - traction[2])
        self.assertTrue(t02 / t01 < 1e-12)

    def test_restart(self):
        # test if restart option works correctly

        # adapt parameters, create solver
        self.parameters['settings']['cores'] = 1
        self.parameters['settings']['flow_iterations'] = 30
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # set displacement
        model_part = interface_input.get_model_part(self.mp_name_in)
        x0, y0, z0 = model_part.x0, model_part.y0, model_part.z0
        dy, dz = self.get_dy_dz(x0, y0, z0)
        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')
        displacement[:, 1] = dy
        displacement[:, 2] = dz

        # run solver for 4 timesteps
        solver.initialize()
        for i in range(4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_1 = solver.get_interface_input()
        interface_y_1 = solver.get_interface_output()
        coords_1 = solver.get_coordinates()[self.mp_name_in]['coords']
        solver.finalize()

        # get data for solver without restart
        interface_output = solver.get_interface_output()
        pressure_1 = interface_output.get_variable_data(self.mp_name_out, 'pressure')
        traction_1 = interface_output.get_variable_data(self.mp_name_out, 'traction')

        # create solver which restarts at timestep 2
        self.parameters['settings']['timestep_start'] = 2
        solver = create_instance(self.parameters)
        interface_input = solver.get_interface_input()

        # run solver for 2 more timesteps
        solver.initialize()
        for i in range(2, 4):
            solver.initialize_solution_step()
            interface_input.set_variable_data(self.mp_name_in, 'displacement', i * displacement)
            solver.solve_solution_step(interface_input)
            solver.finalize_solution_step()
        interface_x_2 = solver.get_interface_input()
        interface_y_2 = solver.get_interface_output()
        coords_2 = solver.get_coordinates()[self.mp_name_in]['coords']
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
        np.testing.assert_allclose(pressure_1, pressure_2, rtol=1e-14)
        np.testing.assert_allclose(traction_1, traction_2, rtol=1e-14)


if __name__ == '__main__':
    unittest.main()
