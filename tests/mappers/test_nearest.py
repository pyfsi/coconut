from coconut.tests.mappers.test_interpolator import Case1D, Case2DCircle, Case3DSphere, Case3DCylinder, Case3DSinc

import unittest


class TestMapperNearest(unittest.TestCase):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.nearest',
                           'settings': {'directions': ['x', 'y', 'z']}}

    def test_mapping_one_dimension(self):
        n_from, n_to = 14, 5
        self.parameters['settings']['directions'] = ['z']

        """
        1D case: square-root grid + linear function

        n_from = 14, n_to = 5 
            => max error = 0.0085
        """
        case = Case1D(n_from, n_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=0.01))
        if self.gui:
            case.plot()

    def test_mapping_circle(self):
        n_from, n_to = 33, 22
        self.parameters['settings']['directions'] = ['x', 'y']

        """
        2D case: circle + linear function

        n_from = 33, n_to = 22 
            => max error = 0.55
        """

        case = Case2DCircle(n_from, n_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=1.))
        if self.gui:
            case.plot()

    def test_mapping_sphere(self):

        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case : sphere + sine function

        n_theta_from, n_phi_from = 40, 20
        n_theta_to, n_phi_to = 50, 21
            => max error = 0.16

        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 0.13
        """

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=0.2))
        if self.gui:
            case.plot()

    def test_mapping_cylinder(self):
        n_x_from, n_theta_from = 18, 20
        n_x_to, n_theta_to = 20, 18
        length = 30.
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case: cylinder + sine function

        n_x_from, n_theta_from = 18, 8
        n_x_to, n_theta_to = 20, 10
        length = 10.
            => max error = 0.14
        """

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=.2))
        if self.gui:
            case.plot()

    def test_mapping_sinc(self):
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case: sinc + linear vector function
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.13
        """

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(self.parameters)
        for tmp in case.check(tolerance=0.2):
            self.assertTrue(tmp)
        if self.gui:
            case.plot()


if __name__ == '__main__':
    unittest.main()
