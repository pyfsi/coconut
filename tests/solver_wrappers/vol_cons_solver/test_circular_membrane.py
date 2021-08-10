from coconut.tools import create_instance, get_solver_env, solver_available

import unittest
import numpy as np
import os
import json


#@unittest.skipUnless(solver_available(f'vol_cons_solver.v1'), f'vol_cons_solver not available')
class TestSolverWrapperVolConsSolverV1(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_circular_membrane', 'parameters.json')

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
        cls.par_solver = parameters['solver_wrappers'][0]

        # if running from this folder
        if os.getcwd() == os.path.realpath(os.path.dirname(__file__)):
            cls.par_solver['settings']['working_directory'] = 'test_circular_membrane'

        cls.folder_path = os.path.join(os.getcwd(), cls.par_solver['settings']['working_directory'])
        cls.delta_t = cls.par_solver['settings']['delta_t']

    @classmethod
    def tearDownClass(cls):
        log_path = os.path.join(cls.folder_path, 'log')
        if os.path.exists(log_path):
            os.remove(log_path)

    def setUp(self):
        self.solver = create_instance(self.par_solver)
        self.mp_name_in = self.solver.settings['interface_input'][0]['model_part']
        self.mp_name_out = self.solver.settings['interface_output'][0]['model_part']
        self.model_part = self.solver.model.get_model_part(self.mp_name_in)

    # noinspection PyMethodMayBeStatic
    def get_displacement(self, points):
        #
        r = np.linalg.norm(points[:,:2], axis=1)
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        # parabolic displacement
        a = 0.1
        dx = np.zeros(points.shape[0])
        dy = np.zeros(points.shape[0])
        dz = (a * np.square(r) - a)
        displacement = np.column_stack((dx, dy, dz))
        return displacement

    # test if nodes are moved to the correct position
    def test_displacement_on_nodes(self):
        x0, y0, z0 = self.model_part.x0, self.model_part.y0, self.model_part.z0
        points = np.stack((x0, y0, z0), axis=1)
        displacement = self.get_displacement(points)
        x = x0 + displacement[:, 0]
        y = y0 + displacement[:, 1]
        z = z0 + displacement[:, 2]
        node_coords_ref = np.stack((x, y, z), axis=1)
        self.solver.get_interface_input().set_variable_data(self.mp_name_in, 'displacement', displacement)

        # update position by iterating once in solver
        self.solver.initialize()
        self.solver.initialize_solution_step()
        self.solver.solve_solution_step(self.solver.get_interface_input())
        self.solver.finalize_solution_step()
        self.solver.finalize()
        node_coords = self.solver.free_surface_updaters[0].get_mesh().points
        np.testing.assert_allclose(node_coords, node_coords_ref, rtol=1e-12)

    #test if same coordinates always gives same pressure & traction
    def test_pressure(self):

        self.solver.initialize()
        self.solver.initialize_solution_step()
        interface_input = self.solver.get_interface_input()

        # set displacement
        x0, y0, z0 = self.model_part.x0, self.model_part.y0, self.model_part.z0
        points = np.stack((x0, y0, z0), axis=1)
        displacement= self.get_displacement(points)
        displacement_list = [displacement, np.zeros_like(displacement), displacement]

        displacement = interface_input.get_variable_data(self.mp_name_in, 'displacement')

        # run solver for three displacements (first one = last one)
        pressure = []
        traction = []
        for disp in displacement_list:
            interface_input.set_variable_data(self.mp_name_in, 'displacement', disp)
            interface_output = self.solver.solve_solution_step(interface_input)
            pressure.append(interface_output.get_variable_data(self.mp_name_out, 'pressure'))
        self.solver.finalize_solution_step()
        self.solver.finalize()

        # check if same position gives same pressure & traction
        np.testing.assert_allclose(pressure[0] , pressure[2] , rtol=1e-12)


        # check if different position gives different pressure & traction
        p01 = np.linalg.norm(pressure[0] - pressure[1])
        p02 = np.linalg.norm(pressure[0] - pressure[2])
        self.assertTrue(p02 / p01 < 1e-12)

if __name__ == '__main__':
    unittest.main()
