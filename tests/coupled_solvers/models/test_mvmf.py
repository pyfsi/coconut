from coconut.tests.coupled_solvers.models import model

import unittest
import numpy as np


class TestModelMVMF(model.TestModel):

    def setUp(self):
        self.q = 2
        self.min_significant = 1.0

        self.parameters = {'type': 'coupled_solvers.models.mvmf',
                           'settings': {
                               'min_significant': self.min_significant,
                               'q': self.q
                           }}
        self.settings = self.parameters['settings']

        super().setUp()

    def test_model(self):
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
        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 1))
        np.testing.assert_array_equal(v[:, 0], self.r2 - self.r1)
        self.assertEqual(w.shape, (self.m, 1))
        np.testing.assert_array_equal(w[:, 0], self.xt2 - self.xt1)

        r.set_interface_data(self.r4)
        xt.set_interface_data(self.xt4)
        self.model.add(r, xt)

        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 2))
        np.testing.assert_array_equal(v[:, 0], self.r4 - self.r2)
        self.assertEqual(w.shape, (self.m, 2))
        np.testing.assert_array_equal(w[:, 0], self.xt4 - self.xt2)

        r.set_interface_data(self.r5)
        xt.set_interface_data(self.xt5)
        self.model.add(r, xt)

        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 3))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack((self.r5 - self.r4, self.r4 - self.r2,
                                                                self.r2 - self.r1)))
        self.assertEqual(w.shape, (self.m, 3))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((self.xt5 - self.xt4, self.xt4 - self.xt2,
                                                                self.xt2 - self.xt1)))

        self.model.filter()

        v = self.model.v
        w = self.model.w
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
        self.model.filter()

        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5, self.r5 - self.r4, self.r4 - self.r2)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack(
            (self.xt8 - self.xt7, self.xt7 - self.xt6, self.xt6 - self.xt5, self.xt5 - self.xt4, self.xt4 - self.xt2)))

        r.set_interface_data(self.r9)
        xt.set_interface_data(self.xt9)
        self.model.add(r, xt)
        self.model.filter()  # adding two modes: filtering at the end is not necessarily the same as filtering each time
        r.set_interface_data(self.r10)
        xt.set_interface_data(self.xt10)
        self.model.add(r, xt)
        self.model.filter()

        v = self.model.v
        w = self.model.w
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
        self.model.filter()

        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 5))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(w.shape, (self.m, 5))
        np.testing.assert_array_equal(w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))

        v1 = self.model.v
        q1, r1 = np.linalg.qr(v1, mode='reduced')
        w1 = self.model.w

        # new solution step
        self.model.finalize_solution_step()
        np.testing.assert_array_equal(self.model.wprev[0], w1)
        np.testing.assert_array_equal(self.model.rrprev[0], r1)
        np.testing.assert_array_equal(self.model.qqprev[0], q1)
        self.model.initialize_solution_step()
        self.assertIsNone(self.model.rref)
        self.assertFalse(self.model.added)
        self.assertEqual(self.model.v.shape, (self.m, 0))
        self.assertEqual(self.model.w.shape, (self.m, 0))
        is_ready = self.model.is_ready()
        self.assertTrue(is_ready)
        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r11)
        xt.set_interface_data(self.xt11)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        self.assertTrue(self.model.added)
        self.assertEqual(self.model.v.shape, (self.m, 0))
        self.assertEqual(self.model.w.shape, (self.m, 0))

        r.set_interface_data(self.r12)
        xt.set_interface_data(self.xt12)
        self.model.add(r, xt)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [2.04216706698691, -8.02212881416246, -4.68760564006761, -1.31217195979005, -8.67687483319990]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        v = self.model.v
        w = self.model.w
        self.assertEqual(v.shape, (self.m, 1))
        np.testing.assert_array_equal(v.T.flatten(), np.hstack(
            (self.r12 - self.r11)))
        self.assertEqual(w.shape, (self.m, 1))
        np.testing.assert_array_equal(w.T.flatten(), np.hstack((self.xt12 - self.xt11)))

        v2 = self.model.v
        q2, r2 = np.linalg.qr(v2, mode='reduced')
        w2 = self.model.w

        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        np.testing.assert_array_equal(self.model.wprev[1], w1)
        np.testing.assert_array_equal(self.model.rrprev[1], r1)
        np.testing.assert_array_equal(self.model.qqprev[1], q1)
        np.testing.assert_array_equal(self.model.wprev[0], w2)
        np.testing.assert_array_equal(self.model.rrprev[0], r2)
        np.testing.assert_array_equal(self.model.qqprev[0], q2)
        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        self.assertEqual(len(self.model.wprev), self.q)
        self.assertEqual(len(self.model.rrprev), self.q)
        self.assertEqual(len(self.model.qqprev), self.q)
        np.testing.assert_array_equal(self.model.wprev[1], w2)
        np.testing.assert_array_equal(self.model.rrprev[1], r2)
        np.testing.assert_array_equal(self.model.qqprev[1], q2)
        np.testing.assert_array_equal(self.model.wprev[0], np.empty((self.m, 0)))
        np.testing.assert_array_equal(self.model.rrprev[0], np.empty((0, 0)))
        np.testing.assert_array_equal(self.model.qqprev[0], np.empty((self.m, 0)))
        self.model.finalize_solution_step()

    def store_old_values(self):
        self.rrprev = self.model.rrprev.copy()
        self.qqprev = self.model.qqprev.copy()
        self.wprev = self.model.wprev.copy()

    def set_new_values(self):
        self.min_significant_new = 0.5
        self.q_new = 1
        self.settings['min_significant'] = self.min_significant_new
        self.settings['q'] = self.q_new

    def check_new_values(self):
        self.assertEqual(self.model.min_significant, self.min_significant_new)
        self.assertEqual(self.model.q, self.q_new)
        np.testing.assert_allclose(np.hstack(self.model.rrprev), np.hstack([self.rrprev[0]]))
        np.testing.assert_allclose(np.hstack(self.model.qqprev), np.hstack([self.qqprev[0]]))
        np.testing.assert_allclose(np.hstack(self.model.wprev), np.hstack([self.wprev[0]]))


if __name__ == '__main__':
    unittest.main()
