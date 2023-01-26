from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import numpy as np
import copy


class TestModel(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        m = cls.m = 5
        dz = 2
        x = 10
        variable = 'area'
        model_part_name = 'wall'
        interface_settings = [{'model_part': model_part_name, 'variables': [variable]}]

        # create model and model_part
        model = data_structure.Model()
        ids = np.arange(0, m)
        x0 = np.zeros(m)
        y0 = np.zeros(m)
        z0 = np.arange(0, m * dz, dz)
        model.create_model_part(model_part_name, x0, y0, z0, ids)

        x_array = np.full((cls.m, 1), x)

        # create interface
        cls.interface = Interface(interface_settings, model)
        cls.interface.set_variable_data(model_part_name, variable, x_array)

        cls.parameters = {}

    def setUp(self):
        self.model = create_instance(self.parameters)
        self.model.size_in = self.model.size_out = self.m
        self.model.out = self.interface.copy()
        self.model.initialize()
        self.model.initialize_solution_step()

        self.r1 = np.array([1, 2, 3, 4, 5])
        self.xt1 = np.array([5, 4, 3, 2, 1])
        self.r2 = np.array([8, 5, 5, 5, 8])
        self.xt2 = np.array([1, 4, 8, 5, 5])
        self.r3 = np.array([7, 5, 6, 4, 3])
        self.xt3 = np.array([5, 4, 2, 5, 1])  # not used everywhere
        self.r4 = np.array([1, 1, 7, 4, 0])
        self.xt4 = np.array([9, 7, 5, 8, 4])
        self.r5 = np.array([9, 5, 10, 6, 4])
        self.xt5 = np.array([5, 1, 2, 3, 9])
        self.r6 = np.array([7, 8, 1, 2, 3])
        self.xt6 = np.array([7, 5, 5, 1, 2])
        self.r7 = np.array([1, 2, 5, 1, 2])
        self.xt7 = np.array([4, 2, 1, 1, 2])
        self.r8 = np.array([6, 3, 9, 0, 3])
        self.xt8 = np.array([3, 1, 2, 3, 9])
        self.r9 = np.array([1, 3, 5, 0, 8])
        self.xt9 = np.array([8, 1, 5, 3, 9])
        self.r10 = np.array([1, 3, -5, 8, 8])
        self.xt10 = np.array([8, -9, 5, 3, -9])
        self.r11 = self.r1
        self.xt11 = self.xt1
        self.r12 = self.r2
        self.xt12 = self.xt2
        self.r13 = self.r10 * 0.95
        self.xt13 = self.xt10

    def test_first_time_step(self):
        r = self.interface.copy()
        xt = self.interface.copy()

        self.assertFalse(self.model.is_ready())

        r.set_interface_data(self.r1)
        xt.set_interface_data(self.xt1)
        self.model.add(r, xt)

        self.assertFalse(self.model.is_ready())

        r.set_interface_data(self.r2)
        xt.set_interface_data(self.xt2)
        self.model.add(r, xt)

        self.assertTrue(self.model.is_ready())

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [4.94444444444445, 0, -6.18055555555556, -3.70833333333333, -4.944444444444452]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r4)
        xt.set_interface_data(self.xt4)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [3.55677154582763, -1.20861833105335, -7.26607387140903, -6.29343365253078, -6.37688098495212]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r5)
        xt.set_interface_data(self.xt5)
        self.model.add(r, xt)
        self.model.filter()

        # test filter_q
        r.set_interface_data(self.r3)
        dr = self.model.filter_q(-1 * r)
        qq, *_ = np.linalg.qr(np.vstack((self.r5 - self.r4, self.r4 - self.r2)).T, mode='reduced')
        dr_in = -self.r3[:np.newaxis]
        dr_sol = dr_in - qq @ (qq.T @ dr_in)
        np.testing.assert_allclose(dr.get_interface_data(), dr_sol)

        r.set_interface_data(self.r6)
        xt.set_interface_data(self.xt6)
        self.model.add(r, xt)
        r.set_interface_data(self.r7)
        xt.set_interface_data(self.xt7)
        self.model.add(r, xt)
        r.set_interface_data(self.r8)
        xt.set_interface_data(self.xt8)
        self.model.add(r, xt)

        # test filter_q
        r.set_interface_data(self.r3)
        dr = self.model.filter_q(-1 * r)
        qq, *_ = np.linalg.qr(np.vstack(
            (self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5, self.r5 - self.r4, self.r4 - self.r2)).T,
            mode='reduced')
        dr_in = -self.r3[:np.newaxis]
        dr_sol = dr_in - qq @ (qq.T @ dr_in)
        np.testing.assert_allclose(dr.get_interface_data(), dr_sol, atol=5e-15)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [7.67847593582868, 21.1203208556145, 21.6136363636358, 23.3649732620315, -1.94652406417126]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r9)
        xt.set_interface_data(self.xt9)
        self.model.add(r, xt)
        self.model.filter()
        r.set_interface_data(self.r10)
        xt.set_interface_data(self.xt10)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r13)
        xt.set_interface_data(self.xt13)
        self.model.add(r, xt)
        r.set_interface_data(self.r3)

        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        # test filter_q
        r.set_interface_data(self.r3)
        dr = self.model.filter_q(-1 * r)
        qq, *_ = np.linalg.qr(np.vstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)).T,
                              mode='reduced')
        dr_in = -self.r3[:np.newaxis]
        dr_sol = dr_in - qq @ (qq.T @ dr_in)
        np.testing.assert_allclose(dr.get_interface_data(), dr_sol, atol=1e-14)

    def test_change_settings_on_restart(self):
        # test if changing settings upon restart works correctly

        r = self.interface.copy()
        xt = self.interface.copy()

        def do_coupling_iteration(r_array, xt_array):
            r.set_interface_data(r_array)
            xt.set_interface_data(xt_array)
            self.model.add(r, xt)
            if self.model.is_ready():
                self.model.predict(r)

        # run for 2 time steps
        for r_array, xt_array in zip((self.r2, self.r3, self.r4, self.r5),
                                     (self.xt2, self.xt3, self.xt4, self.xt5)):
            do_coupling_iteration(r_array, xt_array)

        self.model.finalize_solution_step()
        self.model.initialize_solution_step()

        for r_array, xt_array in zip((self.r6, self.r7, self.r8),
                                     (self.xt6, self.xt7, self.xt8)):
            do_coupling_iteration(r_array, xt_array)

        self.model.finalize_solution_step()
        self.store_old_values()
        self.restart_data = copy.deepcopy(self.model.save_restart_data())
        self.parameters_old = self.parameters.copy()
        self.model.finalize()

        # create restarted model
        self.set_new_values()
        self.model = create_instance(self.parameters)
        self.model.check_restart_data(self.parameters_old)
        self.model.size_in = self.model.size_out = self.m
        self.model.out = self.interface.copy()
        self.model.initialize()
        self.model.restart(self.restart_data)
        self.model.initialize_solution_step()

        self.check_new_values()

        # run for 2 time steps
        r = self.interface.copy()
        xt = self.interface.copy()

        for r_array, xt_array in zip((self.r9, self.r10, self.r11),
                                     (self.xt9, self.xt10, self.xt11)):
            do_coupling_iteration(r_array, xt_array)

        self.model.finalize_solution_step()
        self.model.initialize_solution_step()

        for r_array, xt_array in zip((self.r12, self.r13),
                                     (self.xt12, self.xt13)):
            do_coupling_iteration(r_array, xt_array)

        self.model.finalize_solution_step()
        self.model.finalize()

    def store_old_values(self):
        pass

    def set_new_values(self):
        pass

    def check_new_values(self):
        pass


if __name__ == '__main__':
    unittest.main()
