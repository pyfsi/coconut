from coconut.tests.mappers.test_interpolator import Case1D, Case2D, Case3DSphere, Case3DCylinder, Case3DSinc

import unittest


class TestMapperRadialBasis(unittest.TestCase):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.radial_basis',
                           'settings': {
                               'directions': ['x', 'y', 'z'], 
                               'parallel': True,
                               'include_polynomial': True
                           }}

    def test_mapping_one_dimension(self):
        n_from, n_to = 14, 5
        self.parameters['settings']['directions'] = ['z']

        """
        1D case: square-root grid + linear function
        
        n_from = 14, n_to = 5
        with linear polynomial
            => max error = 5.6e-17
        without linear polynomial
            => max error = 6.1e-7 (sh_par 2: 5.8e-5, sh_par 1000: 1.2e-7)
        """

        case = Case1D(n_from, n_to)

        # with polynomial
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=1e-16))
        if self.gui:
            case.plot(add_text='with polynomial')

        # without polynomial
        self.parameters['settings']['include_polynomial'] = False
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=7e-7))
        if self.gui:
            case.plot(add_text='without polynomial')

    def test_mapping_circle(self):
        n_from, n_to = 33, 22
        self.parameters['settings']['directions'] = ['x', 'y']

        """
        2D case: circle + linear function
        
        n_from = 33, n_to = 22
        with linear polynomial
            => max error = 1.8e-15
        without linear polynomial
            => max error = 1.4e-5 (sh_par 2: 0.003, sh_par 1000: 2.7e-6)
        """

        case = Case2D(n_from, n_to)

        # with polynomial
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=2e-15))
        if self.gui:
            case.plot(add_text='with polynomial')

        # without polynomial
        self.parameters['settings']['include_polynomial'] = False
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=2e-5))
        if self.gui:
            case.plot(add_text='without polynomial')

    def test_mapping_sphere(self):
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        self.parameters['settings']['directions'] = ['x', 'y', 'z']

        """
        3D case: sphere + sine function
        
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 1.0e-4 (sh_par 2: 3.2e-4, sh_par 1000: 1.0e-4 + warning)
        """

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=2e-4))
        if self.gui:
            case.plot()

    def test_mapping_cylinder(self):
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
        self.parameters['settings']['directions'] = ['x', 'y', 'z']
        self.parameters['settings']['scaling'] = [.01, 1, 1]  # bad result without scaling

        """
        3D case: cylinder + sine function
        
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
            => max error = 2.8e-4 (sh_par 2: 4.7-e4, sh_par 1000: 2.8e-4)
        """

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(self.parameters)
        self.assertTrue(case.check(tolerance=3e-4))
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
        with polynomial
            => max error = 6.7e-16
        without polynomial
            => max error = 0.0024 (sh_par 2: 0.05, sh_par 1000: 4.8e-5)
            
        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 7e-6 (sh_par 2: 0.0024, sh_par 1000: 1.4e-6)
        """

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)

        # with polynomial
        case.map(self.parameters)
        for comp in case.check(tolerance=1e-15):
            self.assertTrue(comp)
        if self.gui:
            case.plot(add_text='with polynomial')

        # without polynomial
        self.parameters['settings']['include_polynomial'] = False
        case.map(self.parameters)
        for comp in case.check(tolerance=0.003):
            self.assertTrue(comp)
        if self.gui:
            case.plot(add_text='without polynomial')


if __name__ == '__main__':
    unittest.main()
