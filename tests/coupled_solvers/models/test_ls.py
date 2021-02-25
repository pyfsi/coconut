from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np


class TestModelLS(unittest.TestCase):

    def test_model_ls(self):
        m = 5
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

        x_array = np.full((m, 1), x)

        # create interface
        interface = Interface(interface_settings, model)
        interface.set_variable_data(model_part_name, variable, x_array)

        # read settings
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_ls.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.load(parameter_file)

        # with reuse

        min_significant = settings["setting1"]["settings"]["min_significant"]
        q = settings["setting1"]["settings"]["q"]

        ls = create_instance(settings["setting1"])
        ls.size_in = ls.size_out = m
        ls.out = interface.copy()
        ls.initialize()
        ls.initialize_solution_step()

        r = interface.copy()
        xt = interface.copy()
        r1 = np.array([1, 2, 3, 4, 5])
        xt1 = np.array([5, 4, 3, 2, 1])
        r2 = np.array([8, 5, 5, 5, 8])
        xt2 = np.array([1, 4, 8, 5, 5])
        r3 = np.array([7, 5, 6, 4, 3])
        r4 = np.array([1, 1, 7, 4, 0])
        xt4 = np.array([9, 7, 5, 8, 4])
        r5 = np.array([9, 5, 10, 6, 4])
        xt5 = np.array([5, 1, 2, 3, 9])
        r6 = np.array([7, 8, 1, 2, 3])
        xt6 = np.array([7, 5, 5, 1, 2])
        r7 = np.array([1, 2, 5, 1, 2])
        xt7 = np.array([4, 2, 1, 1, 2])
        r8 = np.array([6, 3, 9, 0, 3])
        xt8 = np.array([3, 1, 2, 3, 9])
        r9 = np.array([1, 3, 5, 0, 8])
        xt9 = np.array([8, 1, 5, 3, 9])
        r10 = np.array([1, 3, -5, 8, 8])
        xt10 = np.array([8, -9, 5, 3, -9])
        r11 = r1
        xt11 = xt1
        r12 = r2
        xt12 = xt2
        r13 = r10 * 0.95
        xt13 = xt10

        is_ready = ls.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(r1)
        xt.set_interface_data(xt1)
        ls.add(r, xt)

        self.assertTrue(ls.added)
        np.testing.assert_array_equal(r1[:,np.newaxis], ls.rref)
        np.testing.assert_array_equal(xt1[:,np.newaxis], ls.xtref)
        is_ready = ls.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(r2)
        xt.set_interface_data(xt2)
        ls.add(r, xt)

        self.assertTrue(ls.added)
        np.testing.assert_array_equal(r2[:,np.newaxis], ls.rref)
        np.testing.assert_array_equal(xt2[:,np.newaxis], ls.xtref)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 1))
        np.testing.assert_array_equal(v[:, 0], r2 - r1)
        self.assertEqual(w.shape, (m, 1))
        np.testing.assert_array_equal(w[:, 0], xt2 - xt1)
        is_ready = ls.is_ready()
        self.assertTrue(is_ready)
        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [4.94444444444445, 0, -6.18055555555556, -3.70833333333333, -4.944444444444452]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(r4)
        xt.set_interface_data(xt4)
        ls.add(r, xt)

        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        np.testing.assert_array_equal(v[:, 0], r4 - r2)
        self.assertEqual(w.shape, (m, 2))
        np.testing.assert_array_equal(w[:, 0], xt4 - xt2)
        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [3.55677154582763, -1.20861833105335, -7.26607387140903, -6.29343365253078, -6.37688098495212]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(r5)
        xt.set_interface_data(xt5)
        ls.add(r, xt)

        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        np.testing.assert_array_equal(v[:, 0], r5 - r4)
        self.assertEqual(w.shape, (m, 2))
        np.testing.assert_array_equal(w[:, 0], xt5 - xt4)
        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [2.17629179331307, 7.43617021276596, 5.80395136778116, 5.96504559270517, -6.89209726443769]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 2))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r5 - r4, r4 - r2)))
        self.assertEqual(w.shape, (m, 2))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt5 - xt4, xt4 - xt2)))

        r.set_interface_data(r6)
        xt.set_interface_data(xt6)
        ls.add(r, xt)
        r.set_interface_data(r7)
        xt.set_interface_data(xt7)
        ls.add(r, xt)
        r.set_interface_data(r8)
        xt.set_interface_data(xt8)
        ls.add(r, xt)

        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [7.67847593582868, 21.1203208556145, 21.6136363636358, 23.3649732620315, -1.94652406417126]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r8 - r7, r7 - r6, r6 - r5, r5 - r4, r4 - r2)))
        self.assertEqual(w.shape, (m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt8 - xt7, xt7 - xt6, xt6 - xt5, xt5 - xt4, xt4 - xt2)))

        r.set_interface_data(r9)
        xt.set_interface_data(xt9)
        ls.add(r, xt)
        r.set_interface_data(r10)
        xt.set_interface_data(xt10)
        ls.add(r, xt)

        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        r.set_interface_data(r13)
        xt.set_interface_data(xt13)
        ls.add(r, xt)
        r.set_interface_data(r3)

        dxt = ls.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        v1 = ls.vcurr
        w1 = ls.wcurr

        # new solution step
        ls.finalize_solution_step()
        np.testing.assert_array_equal(ls.vprev[0].flatten(), v1.flatten())
        np.testing.assert_array_equal(ls.wprev[0].flatten(), w1.flatten())
        ls.initialize_solution_step()
        self.assertIsNone(ls.rref)
        self.assertFalse(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))
        is_ready = ls.is_ready()
        self.assertTrue(is_ready)

        r.set_interface_data(r11)
        xt.set_interface_data(xt11)
        ls.add(r, xt)

        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        self.assertTrue(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))

        r.set_interface_data(r12)
        xt.set_interface_data(xt12)
        ls.add(r, xt)

        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [-8.52029914529913, -3.96866096866096, -0.505163817663819, -0.426103988603991, -14.1239316239316]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r12 - r11, r10 - r9, r9 - r8, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt12 - xt11, xt10 - xt9, xt9 - xt8, xt7 - xt6, xt6 - xt5)))

        v2 = ls.vcurr
        w2 = ls.wcurr

        # new solution step
        ls.finalize_solution_step()
        np.testing.assert_array_equal(np.hstack(ls.vprev).flatten(), np.hstack([v2, v1[:, :2], v1[:, 3:]]).flatten())
        np.testing.assert_array_equal(np.hstack(ls.wprev).flatten(), np.hstack([w2, w1[:, :2], w1[:, 3:]]).flatten())
        ls.initialize_solution_step()

        # new solution step
        ls.finalize_solution_step()
        np.testing.assert_array_equal(np.hstack(ls.vprev).flatten(), np.hstack([np.empty((m, 0)), v2]).flatten())
        np.testing.assert_array_equal(np.hstack(ls.wprev).flatten(), np.hstack([np.empty((m, 0)), w2]).flatten())
        self.assertEqual(len(ls.vprev), q)
        self.assertEqual(len(ls.wprev), q)
        ls.initialize_solution_step()

        # without reuse

        min_significant = settings["setting2"]["settings"]["min_significant"]
        q = settings["setting2"]["settings"]["q"]

        ls = create_instance(settings["setting2"])
        ls.size_in = ls.size_out = m
        ls.out = interface.copy()
        ls.initialize()
        ls.initialize_solution_step()

        r.set_interface_data(r1)
        xt.set_interface_data(xt1)
        ls.add(r, xt)
        r.set_interface_data(r2)
        xt.set_interface_data(xt2)
        ls.add(r, xt)
        r.set_interface_data(r4)
        xt.set_interface_data(xt4)
        ls.add(r, xt)
        r.set_interface_data(r5)
        xt.set_interface_data(xt5)
        ls.add(r, xt)
        r.set_interface_data(r6)
        xt.set_interface_data(xt6)
        ls.add(r, xt)
        r.set_interface_data(r7)
        xt.set_interface_data(xt7)
        ls.add(r, xt)
        r.set_interface_data(r8)
        xt.set_interface_data(xt8)
        ls.add(r, xt)
        r.set_interface_data(r9)
        xt.set_interface_data(xt9)
        ls.add(r, xt)
        r.set_interface_data(r10)
        xt.set_interface_data(xt10)
        ls.add(r, xt)

        r.set_interface_data(r3)
        dxt = ls.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((ls.vcurr, np.hstack(ls.vprev)))
        w = np.hstack((ls.wcurr, np.hstack(ls.wprev)))
        self.assertEqual(v.shape, (m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((r10 - r9, r9 - r8, r8 - r7, r7 - r6, r6 - r5)))
        self.assertEqual(w.shape, (m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((xt10 - xt9, xt9 - xt8, xt8 - xt7, xt7 - xt6, xt6 - xt5)))

        # new solution step
        ls.finalize_solution_step()
        np.testing.assert_array_equal(np.hstack(ls.vprev).flatten(), np.hstack((np.empty((m, 0)))).flatten())
        np.testing.assert_array_equal(np.hstack(ls.wprev).flatten(), np.hstack((np.empty((m, 0)))).flatten())
        ls.initialize_solution_step()
        self.assertIsNone(ls.rref)
        self.assertFalse(ls.added)
        self.assertEqual(ls.vcurr.shape, (m, 0))
        self.assertEqual(ls.wcurr.shape, (m, 0))
        is_ready = ls.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(r11)
        xt.set_interface_data(xt11)
        ls.add(r, xt)

        is_ready = ls.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(r12)
        xt.set_interface_data(xt12)
        ls.add(r, xt)

        is_ready = ls.is_ready()
        self.assertTrue(is_ready)


if __name__ == '__main__':
    unittest.main()
