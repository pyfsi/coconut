from coconut.tests.coupled_solvers.models import model

import unittest
import numpy as np


class TestModelLS(model.TestModel):

    def setUp(self):
        self.parameters['type'] = 'coupled_solvers.models.ls'
        self.q = 2
        self.settings['q'] = self.q
        super().setUp()

    def test_ls_with_reuse(self):
        r = self.interface.copy()
        xt = self.interface.copy()

        r.set_interface_data(self.r1)
        xt.set_interface_data(self.xt1)
        self.model.add(r, xt)

        self.assertTrue(self.model.added)
        np.testing.assert_array_equal(self.r1[:, np.newaxis], self.model.rref)
        np.testing.assert_array_equal(self.xt1[:, np.newaxis], self.model.xtref)

        r.set_interface_data(self.r2)
        xt.set_interface_data(self.xt2)
        self.model.add(r, xt)

        self.assertTrue(self.model.added)
        np.testing.assert_array_equal(self.r2[:, np.newaxis], self.model.rref)
        np.testing.assert_array_equal(self.xt2[:, np.newaxis], self.model.xtref)
        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 1))
        np.testing.assert_array_equal(v[:, 0], self.r2 - self.r1)
        self.assertEqual(w.shape, (self.m, 1))
        np.testing.assert_array_equal(w[:, 0], self.xt2 - self.xt1)

        r.set_interface_data(self.r4)
        xt.set_interface_data(self.xt4)
        self.model.add(r, xt)

        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 2))
        np.testing.assert_array_equal(v[:, 0], self.r4 - self.r2)
        self.assertEqual(w.shape, (self.m, 2))
        np.testing.assert_array_equal(w[:, 0], self.xt4 - self.xt2)

        # test mode limiting
        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r, modes=0)  # use zero modes: return 0
        np.testing.assert_allclose(dxt.get_interface_data(), np.zeros(self.m))
        dxt = self.model.predict(-1 * r, modes=1)  # use one mode: only oldest (first mode)
        dxt_sol = [4.94444444444445, 0, -6.18055555555556, -3.70833333333333, -4.944444444444452]
        # dxt_sol = [5.07462687, 1.90298507, -1.90298507, 1.90298507, -0.63432836]  # only most recent
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r5)
        xt.set_interface_data(self.xt5)
        self.model.add(r, xt)

        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 2))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((self.r5 - self.r4, self.r4 - self.r2)))
        self.assertEqual(w.shape, (self.m, 2))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((self.xt5 - self.xt4, self.xt4 - self.xt2)))

        r.set_interface_data(self.r6)
        xt.set_interface_data(self.xt6)
        self.model.add(r, xt)
        r.set_interface_data(self.r7)
        xt.set_interface_data(self.xt7)
        self.model.add(r, xt)
        r.set_interface_data(self.r8)
        xt.set_interface_data(self.xt8)
        self.model.add(r, xt)

        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5, self.r5 - self.r4, self.r4 - self.r2)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack(
            (self.xt8 - self.xt7, self.xt7 - self.xt6, self.xt6 - self.xt5, self.xt5 - self.xt4, self.xt4 - self.xt2)))

        r.set_interface_data(self.r9)
        xt.set_interface_data(self.xt9)
        self.model.add(r, xt)
        r.set_interface_data(self.r10)
        xt.set_interface_data(self.xt10)
        self.model.add(r, xt)

        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))

        r.set_interface_data(self.r13)
        xt.set_interface_data(self.xt13)
        self.model.add(r, xt)

        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))

        v1 = self.model.vcurr
        w1 = self.model.wcurr

        # new solution step
        self.model.finalize_solution_step()
        np.testing.assert_array_equal(self.model.vprev[0].flatten(), v1.flatten())
        np.testing.assert_array_equal(self.model.wprev[0].flatten(), w1.flatten())
        self.model.initialize_solution_step()
        self.assertIsNone(self.model.rref)
        self.assertFalse(self.model.added)
        self.assertEqual(self.model.vcurr.shape, (self.m, 0))
        self.assertEqual(self.model.wcurr.shape, (self.m, 0))
        is_ready = self.model.is_ready()
        self.assertTrue(is_ready)

        r.set_interface_data(self.r11)
        xt.set_interface_data(self.xt11)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        self.assertTrue(self.model.added)
        self.assertEqual(self.model.vcurr.shape, (self.m, 0))
        self.assertEqual(self.model.wcurr.shape, (self.m, 0))

        r.set_interface_data(self.r12)
        xt.set_interface_data(self.xt12)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-8.52029914529913, -3.96866096866096, -0.505163817663819, -0.426103988603991, -14.1239316239316]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r12 - self.r11, self.r10 - self.r9, self.r9 - self.r8, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(),
                                      np.hstack((self.xt12 - self.xt11, self.xt10 - self.xt9, self.xt9 - self.xt8,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))

        v2 = self.model.vcurr
        w2 = self.model.wcurr
        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        np.testing.assert_array_equal(np.hstack(self.model.vprev).flatten(),
                                      np.hstack([v2, v1[:, :2], v1[:, 3:]]).flatten())
        np.testing.assert_array_equal(np.hstack(self.model.wprev).flatten(),
                                      np.hstack([w2, w1[:, :2], w1[:, 3:]]).flatten())
        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        np.testing.assert_array_equal(np.hstack(self.model.vprev).flatten(),
                                      np.hstack([np.empty((self.m, 0)), v2]).flatten())
        np.testing.assert_array_equal(np.hstack(self.model.wprev).flatten(),
                                      np.hstack([np.empty((self.m, 0)), w2]).flatten())
        self.assertEqual(len(self.model.vprev), self.q)
        self.assertEqual(len(self.model.wprev), self.q)

    def test_ls_without_reuse(self):
        self.parameters['type'] = 'coupled_solvers.models.ls'
        self.q = 0
        self.settings['q'] = self.q
        super().setUp()

        r = self.interface.copy()
        xt = self.interface.copy()

        r.set_interface_data(self.r1)
        xt.set_interface_data(self.xt1)
        self.model.add(r, xt)
        r.set_interface_data(self.r2)
        xt.set_interface_data(self.xt2)
        self.model.add(r, xt)
        r.set_interface_data(self.r4)
        xt.set_interface_data(self.xt4)
        self.model.add(r, xt)
        r.set_interface_data(self.r5)
        xt.set_interface_data(self.xt5)
        self.model.add(r, xt)
        r.set_interface_data(self.r6)
        xt.set_interface_data(self.xt6)
        self.model.add(r, xt)
        r.set_interface_data(self.r7)
        xt.set_interface_data(self.xt7)
        self.model.add(r, xt)
        r.set_interface_data(self.r8)
        xt.set_interface_data(self.xt8)
        self.model.add(r, xt)
        r.set_interface_data(self.r9)
        xt.set_interface_data(self.xt9)
        self.model.add(r, xt)
        r.set_interface_data(self.r10)
        xt.set_interface_data(self.xt10)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = np.hstack((self.model.vcurr, np.hstack(self.model.vprev)))
        w = np.hstack((self.model.wcurr, np.hstack(self.model.wprev)))
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))
        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        np.testing.assert_array_equal(np.hstack(self.model.vprev).flatten(),
                                      np.hstack((np.empty((self.m, 0)))).flatten())
        np.testing.assert_array_equal(np.hstack(self.model.wprev).flatten(),
                                      np.hstack((np.empty((self.m, 0)))).flatten())
        self.assertIsNone(self.model.rref)
        self.assertFalse(self.model.added)
        self.assertEqual(self.model.vcurr.shape, (self.m, 0))
        self.assertEqual(self.model.wcurr.shape, (self.m, 0))
        is_ready = self.model.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(self.r11)
        xt.set_interface_data(self.xt11)
        self.model.add(r, xt)

        is_ready = self.model.is_ready()
        self.assertFalse(is_ready)

        r.set_interface_data(self.r12)
        xt.set_interface_data(self.xt12)
        self.model.add(r, xt)

        is_ready = self.model.is_ready()
        self.assertTrue(is_ready)


if __name__ == '__main__':
    unittest.main()
