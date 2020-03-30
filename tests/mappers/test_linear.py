from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.tests.mappers.test_nearest import Case1D, Case2D, Case3DSphere, Case3DSinc
from coconut.coupling_components.mappers.linear import *

import os
import numpy as np


class TestMapperLinear(KratosUnittest.TestCase):
    def test_mapper_linear(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_linear.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_mapper = parameters['mapper']

        gui = 0  # *** gui gives problems when running all tests?

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 0
        """
        n_from, n_to = 14, 5
        par_mapper['settings'].SetArray('directions', ['Z'])

        case = Case1D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=1e-12))
        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 0.032
        """
        n_from, n_to = 33, 22
        par_mapper['settings'].SetArray('directions', ['X', 'Y'])

        case = Case2D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=0.05))
        if gui:
            case.plot()

        # 3D case: sphere + sine function
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 0.016
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=0.03))
        if gui:
            case.plot()

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.13

        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 0.082
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(par_mapper)
        for tmp in case.check(tolerance=0.2):
            self.assertTrue(tmp)
        if gui:
            case.plot()


        def assertArrayAlmostEqual(a, b, delta):
            for el_a, el_b in zip(a.flatten(), b.flatten()):
                self.assertAlmostEqual(el_a, el_b, delta=delta)

        # test function line_interpolation_coeff
        if True:
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

        # test function get_coeffs_1d_2d
        if True:
            # 1D: linear, nearest neighbour
            coords_from = np.array([[1.], [0.]])

            coord_to = np.array([0.6])
            coeffs = get_coeffs_1d_2d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[.6, .4]]), 1e-15)

            coord_to = np.array([1.4])
            coeffs = get_coeffs_1d_2d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[1., 0.]]), 1e-15)

            # 2D: linear, nearest neighbour
            coords_from = np.array([[1., 1.],
                                    [0., 0.]])

            coord_to = np.array([0., 1.])
            coeffs = get_coeffs_1d_2d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[.5, .5]]), 1e-15)

            coord_to = np.array([2., 1.])
            coeffs = get_coeffs_1d_2d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[1., 0.]]), 1e-15)

        # test function triangle_area
        if True:
            P_1 = np.array([0., 0., 0.])
            P_2 = np.array([1., 0., 0.])
            P_3 = np.array([1., 1., 0.])
            self.assertAlmostEqual(triangle_area(P_1, P_2, P_3), .5, delta=1e-15)

        # test function degenerate_triangle
        if True:
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

        # test function project_on_triangle
        if True:
            P_1 = np.array([0., 0., 0.])
            P_2 = np.array([1., 0., 0.])
            P_3 = np.array([1., 1., 0.])

            P_0 = np.array([.5, .5, 1.])
            P_p = project_on_triangle(P_0, P_1, P_2, P_3)
            assertArrayAlmostEqual(P_p, np.array([[.5, .5, 0.]]), 1e-15)

            P_0 = np.array([0., 1., 1.])
            P_p = project_on_triangle(P_0, P_1, P_2, P_3)
            assertArrayAlmostEqual(P_p, np.array([[0., 1., 0.]]), 1e-15)

        # test function point_on_triangle
        if True:
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

        # test function get_coeffs_3d
        if True:
            coords_from = np.array([[0., 0., 0.],
                                    [1., 0., 0.],
                                    [1., 1., 0.]])

            coord_to = np.array([.75, .25, 0.])
            coeffs = get_coeffs_3d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[.25, .5, .25]]), 1e-15)

            coord_to = np.array([.5, .25, 0.])
            coeffs = get_coeffs_3d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[.5, .25, .25]]), 1e-15)

            coord_to = np.array([.5, -.1, 0.])
            coeffs = get_coeffs_3d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[.5, .5, 0.]]), 1e-15)

            coord_to = np.array([-.1, -.1, 0.])
            coeffs = get_coeffs_3d(coords_from, coord_to)
            assertArrayAlmostEqual(coeffs, np.array([[1., 0., 0.]]), 1e-15)


if __name__ == '__main__':
    KratosUnittest.main()