from coconut.tests.mappers.test_interpolator import Case1D, Case2D, Case3DSphere, Case3DCylinder, Case3DSinc

import unittest
import json


class TestMapperRadialBasis(unittest.TestCase):

    def test_mapper_radial_basis(self):
        with open('mappers/test_radial_basis.json') as parameter_file:
            parameters = json.load(parameter_file)
        par_mapper = parameters['mapper']

        gui = 0

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 6.1e-7 (sh_par 2: 5.8e-5, sh_par 1000: 1.2e-7)
        """
        n_from, n_to = 14, 5
        par_mapper['settings']['directions'] = ['z']

        case = Case1D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=7e-7))
        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 1.4e-5 (sh_par 2: 0.003, sh_par 1000: 2.7e-6)
        """
        n_from, n_to = 33, 22
        par_mapper['settings']['directions'] = ['x', 'y']

        case = Case2D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=2e-5))
        if gui:
            case.plot()

        # 3D case: sphere + sine function
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 1.0e-4 (sh_par 2: 3.2e-4, sh_par 1000: 1.0e-4 + warning)
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        par_mapper['settings']['directions'] = ['x', 'y', 'z']

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=2e-4))
        if gui:
            case.plot()

        # 3D case: cylinder + sine function
        """
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
            => max error = 2.8e-4 (sh_par 2: 4.7-e4, sh_par 1000: 2.8e-4)
        """
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
        par_mapper['settings']['directions'] = ['x', 'y', 'z']
        par_mapper['settings']['scaling'] = [.01, 1, 1]  # bad result without scaling

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=3e-4))
        if gui:
            case.plot()
        par_mapper['settings'].pop('scaling')

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.0024 (sh_par 2: 0.05, sh_par 1000: 4.8e-5)
            
        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 7e-6 (sh_par 2: 0.0024, sh_par 1000: 1.4e-6)
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        par_mapper['settings']['directions'] = ['x', 'y', 'z']

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(par_mapper)
        for tmp in case.check(tolerance=0.003):
            self.assertTrue(tmp)
        if gui:
            case.plot()


if __name__ == '__main__':
    unittest.main()
