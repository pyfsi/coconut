from coconut.tests.coupled_solvers.models import model

import unittest
import numpy as np


class TestModelMV(model.TestModel):

    def setUp(self):
        self.min_significant = 1.0

        self.parameters = {'type': 'coupled_solvers.models.mv',
                           'settings': {
                               'min_significant': self.min_significant
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
        self.assertIsNone(self.model.ncurr)

        r.set_interface_data(self.r2)
        xt.set_interface_data(self.xt2)
        self.model.add(r, xt)

        self.assertTrue(self.model.added)
        np.testing.assert_array_equal(self.r2[:, np.newaxis], self.model.rref)
        np.testing.assert_array_equal(self.xt2[:, np.newaxis], self.model.xtref)
        self.assertEqual(self.model.v.shape, (self.m, 1))
        np.testing.assert_array_equal(self.model.v[:, 0], self.r2 - self.r1)
        self.assertEqual(self.model.w.shape, (self.m, 1))
        np.testing.assert_array_equal(self.model.w[:, 0], self.xt2 - self.xt1)
        n_sol = [-0.388888888888889, -0.166666666666667, -0.111111111111111, -0.0555555555555556, -0.166666666666667,
                 0, 0, 0, 0, 0,
                 0.486111111111111, 0.208333333333333, 0.138888888888889, 0.0694444444444445, 0.208333333333333,
                 0.291666666666667, 0.125000000000000, 0.0833333333333333, 0.0416666666666667, 0.125000000000000,
                 0.388888888888889, 0.166666666666667, 0.111111111111111, 0.0555555555555556, 0.166666666666667]
        np.testing.assert_allclose(self.model.ncurr.flatten(), n_sol)
        np.testing.assert_array_equal(self.model.nprev.flatten(), np.zeros((self.m, self.m)).flatten())

        r.set_interface_data(self.r4)
        xt.set_interface_data(self.xt4)
        self.model.add(r, xt)

        self.assertEqual(self.model.v.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.v[:, 0], self.r4 - self.r2)
        self.assertEqual(self.model.w.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.w[:, 0], self.xt4 - self.xt2)
        n_sol = [-0.306429548563612, -0.216142270861833, 0.251709986320109, -0.0437756497948016, -0.555403556771546,
                 0.0718194254445965, -0.0430916552667578, 0.316005471956224, 0.0102599179206566, -0.338577291381669,
                 0.550615595075240, 0.169630642954856, 0.422708618331053, 0.0786593707250342, -0.0957592339261285,
                 0.445280437756498, 0.0328317373461012, 0.759233926128591, 0.0636114911080711, -0.599179206566348,
                 0.474008207934337, 0.115595075239398, 0.485636114911081, 0.0677154582763338, -0.234610123119015]
        np.testing.assert_allclose(self.model.ncurr.flatten(), n_sol)

        r.set_interface_data(self.r5)
        xt.set_interface_data(self.xt5)
        self.model.add(r, xt)

        self.assertEqual(self.model.v.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.v[:, 0], self.r5 - self.r4)
        self.assertEqual(self.model.w.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.w[:, 0], self.xt5 - self.xt4)
        n_sol = [-0.258792878853669, -0.180633955709944, 0.376899696048632, 0.0121580547112462, -0.590534085974816,
                 -0.460486322188450, -0.200607902735562, -0.446808510638298, -0.159574468085106, 0.0364741641337384,
                 -0.266391663048198, -0.0651324359531047, -0.729483282674772, -0.168693009118541, 0.479374728614850,
                 -0.379722101606600, -0.171081198436822, -0.316109422492401, -0.123100303951368, -0.0208423795049935,
                 0.395788102475033, 0.155449413808076, 0.541033434650456, 0.162613981762918, -0.184107685627443]
        np.testing.assert_allclose(self.model.ncurr.flatten(), n_sol)
        self.assertEqual(self.model.v.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.v.T.flatten(), np.hstack((self.r5 - self.r4, self.r4 - self.r2)))
        self.assertEqual(self.model.w.shape, (self.m, 2))
        np.testing.assert_array_equal(self.model.w.T.flatten(), np.hstack((self.xt5 - self.xt4, self.xt4 - self.xt2)))

        r.set_interface_data(self.r6)
        xt.set_interface_data(self.xt6)
        self.model.add(r, xt)
        r.set_interface_data(self.r7)
        xt.set_interface_data(self.xt7)
        self.model.add(r, xt)
        r.set_interface_data(self.r8)
        xt.set_interface_data(self.xt8)
        self.model.add(r, xt)

        n_sol = [1.59692513368984, -1.67045454545455, -1.19117647058823, 0.612967914438510, -1.93649732620321,
                 3.87433155080215, -5.22727272727274, -3.11764705882354, 0.660427807486646, -2.01336898395721,
                 4.65909090909091, -6.15909090909093, -3.50000000000001, 0.568181818181834, -1.56818181818181,
                 5.05213903743316, -7.02272727272729, -3.32352941176471, 0.736631016042804, -2.20721925133689,
                 1.72192513368984, -1.54545454545455, 0.0588235294117627, -0.262032085561494, -0.561497326203206]
        np.testing.assert_allclose(self.model.ncurr.flatten(), n_sol)
        self.assertEqual(self.model.v.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.v.T.flatten(), np.hstack(
            (self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5, self.r5 - self.r4, self.r4 - self.r2)))
        self.assertEqual(self.model.w.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.w.T.flatten(),
                                      np.hstack((self.xt8 - self.xt7, self.xt7 - self.xt6, self.xt6 - self.xt5,
                                                 self.xt5 - self.xt4, self.xt4 - self.xt2)))

        r.set_interface_data(self.r9)
        xt.set_interface_data(self.xt9)
        self.model.add(r, xt)
        r.set_interface_data(self.r10)
        xt.set_interface_data(self.xt10)
        self.model.add(r, xt)

        self.assertEqual(self.model.v.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(self.model.w.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))

        r.set_interface_data(self.r13)
        xt.set_interface_data(self.xt13)
        self.model.add(r, xt)

        self.assertEqual(self.model.v.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.v.T.flatten(), np.hstack(
            (self.r10 - self.r9, self.r9 - self.r8, self.r8 - self.r7, self.r7 - self.r6, self.r6 - self.r5)))
        self.assertEqual(self.model.w.shape, (self.m, 5))
        np.testing.assert_array_equal(self.model.w.T.flatten(),
                                      np.hstack((self.xt10 - self.xt9, self.xt9 - self.xt8, self.xt8 - self.xt7,
                                                 self.xt7 - self.xt6, self.xt6 - self.xt5)))
        v = self.model.v
        w = self.model.w

        nprev = w @ np.linalg.inv(v.T @ v) @ v.T
        self.model.finalize_solution_step()

        # new solution step
        self.model.initialize_solution_step()
        np.testing.assert_array_equal(self.model.nprev.flatten(), nprev.flatten())
        self.assertIsNone(self.model.rref)
        self.assertFalse(self.model.added)
        self.assertEqual(self.model.v.shape, (self.m, 0))
        self.assertEqual(self.model.w.shape, (self.m, 0))
        is_ready = self.model.is_ready()
        self.assertTrue(is_ready)
        np.testing.assert_allclose(self.model.ncurr.flatten(), nprev.flatten())
        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [-4.19875900720576, -5.62710168134507, -2.21637309847878, -0.788630904723781, -11.8953162530024]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)

        r.set_interface_data(self.r11)
        xt.set_interface_data(self.xt11)
        self.model.add(r, xt)

        self.assertTrue(self.model.added)
        self.assertEqual(self.model.v.shape, (self.m, 0))
        self.assertEqual(self.model.w.shape, (self.m, 0))

        r.set_interface_data(self.r12)
        xt.set_interface_data(self.xt12)
        self.model.add(r, xt)

        n_sol = [-1.07953029089939, 0.852682145716574, -0.00813984520949948, 0.0950093408059783, 0.306645316253002,
                 -1.08484565430122, 2.54250066720043, 0.878480562227564, -0.192264923049552, -0.532759540966108,
                 0.315863802152833, 0.444989324793168, -0.139022328974290, -0.215427897873855, 0.649152655457700,
                 0.519159772262255, -0.566019482252470, -0.0682990837114143, -0.0941975802864517, 0.431578596210302,
                 -0.394248732319190, 0.855284227381908, 1.23271950894049, -0.654857219108619, 0.794435548438751]
        np.testing.assert_allclose(self.model.ncurr.flatten(), n_sol)

        r.set_interface_data(self.r3)
        dxt = self.model.predict(-1 * r)
        dxt_sol = [2.04216706698691, -8.02212881416246, -4.68760564006761, -1.31217195979005, -8.67687483319990]
        np.testing.assert_allclose(dxt.get_interface_data(), dxt_sol)
        self.assertEqual(self.model.v.shape, (self.m, 1))
        np.testing.assert_array_equal(self.model.v, self.r12[:, np.newaxis] - self.r11[:, np.newaxis])
        self.assertEqual(self.model.w.shape, (self.m, 1))
        np.testing.assert_array_equal(self.model.w, self.xt12[:, np.newaxis] - self.xt11[:, np.newaxis])

    def store_old_values(self):
        self.nprev = self.model.nprev.copy()

    def set_new_values(self):
        self.min_significant_new = 0.5
        self.settings['min_significant'] = self.min_significant_new

    def check_new_values(self):
        self.assertEqual(self.model.min_significant, self.min_significant_new)
        np.testing.assert_allclose(self.model.nprev, self.nprev)


if __name__ == '__main__':
    unittest.main()
