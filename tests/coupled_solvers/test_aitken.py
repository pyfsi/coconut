from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance, cd
from coconut.tests.coupled_solvers import coupled_solver

import unittest
import numpy as np


class TestCoupledSolverAitken(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_aitken.json'

    def test_coupled_solver_aitken(self):
        m = 10
        dz = 2
        r = 0.1
        x = 10
        xt0 = 10.5
        xt1 = 10.2
        xt2 = 10.1
        xt3 = 10.7
        xt4 = 9.9
        variable = 'area'  # TODO: does not match JSON setings
        model_part_name = 'wall'
        interface_settings = [{'model_part': model_part_name, 'variables': [variable]}]

        # create model and model_part
        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)
        model.create_model_part(model_part_name, x0, y0, z0, ids)
        
        # create interface
        interface = Interface(interface_settings, model)

        omega_max = self.parameters['settings']['omega_max']

        with cd(self.working_dir):
            coupled_solver = create_instance(self.parameters)
            coupled_solver.initialize()
        coupled_solver.initialize_solution_step()

        interface_r = interface.copy()
        interface_r.set_interface_data(np.full(m, r))
        interface_x = interface.copy()
        interface_x.set_interface_data(x * np.ones(m))
        interface_xt0 = interface.copy()
        interface_xt0.set_interface_data(xt0 * np.ones(m))
        interface_xt1 = interface.copy()
        interface_xt1.set_interface_data(xt1 * np.ones(m))
        interface_xt2 = interface.copy()
        interface_xt2.set_interface_data(xt2 * np.ones(m))
        interface_xt3 = interface.copy()
        interface_xt3.set_interface_data(xt3 * np.ones(m))
        interface_xt4 = interface.copy()
        interface_xt4.set_interface_data(xt4 * np.ones(m))

        # test value of self.added
        is_ready = coupled_solver.is_ready()
        self.assertFalse(is_ready)

        # test update()
        coupled_solver.update(interface_x, interface_xt0)
        is_ready = coupled_solver.is_ready()
        self.assertTrue(is_ready)
        omega = coupled_solver.omega
        self.assertEqual(omega, omega_max)
        coupled_solver.update(interface_x, interface_xt1)
        is_ready = coupled_solver.is_ready()
        self.assertTrue(is_ready)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, 5 / 3 * omega_max, 10)
        coupled_solver.update(interface_x, interface_xt2)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, 10 / 3 * omega_max, 10)

        # test predict()
        interface_dx = coupled_solver.predict(interface_r)
        omega = interface_dx.get_interface_data()[0] / r
        self.assertEqual(omega, coupled_solver.omega)
        coupled_solver.predict(interface_r)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, 10 / 3 * omega_max, 10)

        # new solution step
        coupled_solver.finalize_solution_step()
        coupled_solver.initialize_solution_step()

        # test value of self.added
        is_ready = coupled_solver.is_ready()
        self.assertFalse(is_ready)

        # test update()
        coupled_solver.update(interface_x, interface_xt0)
        is_ready = coupled_solver.is_ready()
        self.assertTrue(is_ready)
        omega = coupled_solver.omega
        self.assertEqual(omega, omega_max)
        coupled_solver.update(interface_x, interface_xt3)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, -5 / 2 * omega_max, 10)

        # new solution step
        coupled_solver.finalize_solution_step()
        coupled_solver.initialize_solution_step()

        # test update()
        coupled_solver.update(interface_x, interface_xt0)
        omega = coupled_solver.omega
        self.assertEqual(omega, -omega_max)
        coupled_solver.update(interface_x, interface_xt4)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, -5 / 6 * omega_max, 10)

        # new solution step
        coupled_solver.finalize_solution_step()
        coupled_solver.initialize_solution_step()

        # test update()
        coupled_solver.update(interface_x, interface_xt0)
        omega = coupled_solver.omega
        self.assertAlmostEqual(omega, -5 / 6 * omega_max, 10)

    def store_old_values(self, coupled_solver):
        self.omega = coupled_solver.omega

    def set_new_values(self):
        self.omega_max_new = 0.6
        self.settings['omega_max'] = self.omega_max_new

    def check_new_values(self, coupled_solver):
        self.assertEqual(coupled_solver.omega, self.omega)
        self.assertEqual(coupled_solver.omega_max, self.omega_max_new)


if __name__ == '__main__':
    unittest.main()
