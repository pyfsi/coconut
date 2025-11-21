from coconut.tests.mappers.test_interpolator import Case1D, Case2DCircle, Case3DSphere, Case3DCylinder, Case3DSinc
from coconut.coupling_components.mappers.linear import *

import unittest
import numpy as np


class TestMapperLinear(unittest.TestCase):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.linear',
                           'settings': {
                               'directions': ['x', 'y', 'z'],
                               'parallel': True
                           }}

    def test_mapping_one_dimension(self):
        n_from, n_to = 14, 5
        self.parameters['settings']['directions'] = ['z']

        """
        1D case: square-root grid + linear function

        n_from = 14, n_to = 5 
            => max error = 0
        """

        case = Case1D(n_from, n_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=1e-12))
        if self.gui:
            case.plot()

    def test_mapping_circle(self):
        n_from, n_to = 33, 22
        self.parameters['settings']['directions'] = ['x', 'y']

        """
        2D case: circle + linear function

        n_from = 33, n_to = 22 
            => max error = 0.032
        """

        case = Case2DCircle(n_from, n_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=0.05))
        if self.gui:
            case.plot()

    def test_mapping_sphere(self):
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case: sphere + sine function

        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 0.016
        """
        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=0.03))
        if self.gui:
            case.plot()

    def test_mapping_cylinder(self):
        n_x_from, n_theta_from = 28, 56
        n_x_to, n_theta_to = 35, 60
        length = 10.
        self.parameters['settings']['directions'] = ['x', 'y', 'z']
        self.parameters['settings']['scaling'] = [.1, 1, 1]  # bad result without scaling

        """
        3D case: cylinder + sine function

        n_x_from, n_theta_from = 28, 56
        n_x_to, n_theta_to = 35, 60
            length = 10.
            => max error = 0.021
        """

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=0.03))
        if self.gui:
            case.plot()
        self.parameters['settings'].pop('scaling')

    def test_mapping_sinc(self):
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case: sinc + linear vector function

        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.13

        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 0.082
        """

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(self.parameters)
        for tmp in case.check(tolerance=0.2):
            self.assertTrue(tmp)
        if self.gui:
            case.plot()

    def test_line_interpolation_coeff(self):
        # 1D: linear, nearest neighbour
        P_1 = np.array([1])
        P_2 = np.array([0])

        P_0 = np.array([0.6])
        c = line_interpolation_coeff(P_0, P_1, P_2)
        self.assertAlmostEqual(c, 0.6, delta=1e-15)

        P_0 = np.array([1.4])
        c = line_interpolation_coeff(P_0, P_1, P_2)
        self.assertAlmostEqual(c, 1., delta=1e-15)

        # 2D: linear, nearest neighbour
        P_1 = np.array([1., 1.])
        P_2 = np.array([0., 0.])

        P_0 = np.array([0., 1.])
        c = line_interpolation_coeff(P_0, P_1, P_2)
        self.assertAlmostEqual(c, 0.5, delta=1e-15)

        P_0 = np.array([2., 1.])
        c = line_interpolation_coeff(P_0, P_1, P_2)
        self.assertAlmostEqual(c, 1., delta=1e-15)

    def test_get_coeffs_1d_2d(self):
        # 1D: linear, nearest neighbour
        coords_from = np.array([[1.], [0.]])

        coord_to = np.array([0.6])
        coeffs = get_coeffs_1d_2d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[.6, .4]]), rtol=1e-15)

        coord_to = np.array([1.4])
        coeffs = get_coeffs_1d_2d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[1., 0.]]), rtol=1e-15)

        # 2D: linear, nearest neighbour
        coords_from = np.array([[1., 1.], [0., 0.]])

        coord_to = np.array([0., 1.])
        coeffs = get_coeffs_1d_2d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[.5, .5]]), rtol=1e-15)

        coord_to = np.array([2., 1.])
        coeffs = get_coeffs_1d_2d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[1., 0.]]), rtol=1e-15)

    def test_triangle_area(self):
        P_1 = np.array([0., 0., 0.])
        P_2 = np.array([1., 0., 0.])
        P_3 = np.array([1., 1., 0.])
        self.assertAlmostEqual(triangle_area(P_1, P_2, P_3), .5, delta=1e-15)

    def test_degenerate_triangle(self):
        P_1 = np.array([0., 0., 0.])
        P_2 = np.array([1., 0., 0.])

        P_3 = np.array([1., 1., 0.])
        self.assertFalse(degenerate_triangle(P_1, P_2, P_3))

        P_3 = np.array([.5, 0., 0.])
        self.assertTrue(degenerate_triangle(P_1, P_2, P_3))

        P_3 = np.array([.5, .01, 0.])
        self.assertTrue(degenerate_triangle(P_1, P_2, P_3))

        P_3 = np.array([0., .01, 0.])
        self.assertTrue(degenerate_triangle(P_1, P_2, P_3))

    def test_project_on_triangle(self):
        P_1 = np.array([0., 0., 0.])
        P_2 = np.array([1., 0., 0.])
        P_3 = np.array([1., 1., 0.])

        P_0 = np.array([.5, .5, 1.])
        P_p = project_on_triangle(P_0, P_1, P_2, P_3)
        np.testing.assert_allclose(P_p, np.array([.5, .5, 0.]), rtol=1e-15)

        P_0 = np.array([0., 1., 1.])
        P_p = project_on_triangle(P_0, P_1, P_2, P_3)
        np.testing.assert_allclose(P_p, np.array([0., 1., 0.]), rtol=1e-15)

    def test_point_on_triangle(self):
        P_1 = np.array([0., 0., 0.])
        P_2 = np.array([1., 0., 0.])
        P_3 = np.array([1., 1., 0.])

        P_p = np.array([.5, .5, 0.])
        self.assertTrue(point_on_triangle(P_p, P_1, P_2, P_3))

        P_p = np.array([0., 1., 0.])
        self.assertFalse(point_on_triangle(P_p, P_1, P_2, P_3))

        P_p = np.array([.5, 1e-6, 0.])
        self.assertTrue(point_on_triangle(P_p, P_1, P_2, P_3))

        P_p = np.array([.5, -1e-6, 0.])
        self.assertFalse(point_on_triangle(P_p, P_1, P_2, P_3))

    def test_get_coeffs_3d(self):
        coords_from = np.array([[0., 0., 0.],
                                [1., 0., 0.],
                                [1., 1., 0.]])

        coord_to = np.array([.75, .25, 0.])
        coeffs = get_coeffs_3d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[.25, .5, .25]]), rtol=1e-15)

        coord_to = np.array([.5, .25, 0.])
        coeffs = get_coeffs_3d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[.5, .25, .25]]), rtol=1e-15)

        coord_to = np.array([.5, -.1, 0.])
        coeffs = get_coeffs_3d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[.5, .5, 0.]]), rtol=1e-15)

        coord_to = np.array([-.1, -.1, 0.])
        coeffs = get_coeffs_3d(coords_from, coord_to)
        np.testing.assert_allclose(coeffs, np.array([[1., 0., 0.]]), rtol=1e-15)


if __name__ == '__main__':
    unittest.main()
