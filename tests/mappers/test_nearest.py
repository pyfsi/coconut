from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.tests.mappers.test_interpolator import Case1D, Case2D, Case3DSphere, Case3DCylinder, Case3DSinc

import os

variables = vars(data_structure)


class TestMapperNearest(KratosUnittest.TestCase):
    def test_mapper_nearest(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_nearest.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_mapper = parameters['mapper']

        gui = 0  # *** gui gives problems when running all tests?

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 0.0085
        """
        n_from, n_to = 14, 5
        par_mapper['settings'].SetArray('directions', ['Z'])

        case = Case1D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=0.01))
        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 0.55
        """
        n_from, n_to = 33, 22
        par_mapper['settings'].SetArray('directions', ['X', 'Y'])

        case = Case2D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=1.))
        if gui:
            case.plot()

        # 3D case : sphere + sine function
        """
        n_theta_from, n_phi_from = 40, 20
        n_theta_to, n_phi_to = 50, 21
            => max error = 0.16

        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 0.13
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=0.2))
        if gui:
            case.plot()

        # 3D case: cylinder + sine function
        """
        n_x_from, n_theta_from = 18, 8
        n_x_to, n_theta_to = 20, 10
        length = 10.
            => max error = 0.14
        """
        n_x_from, n_theta_from = 18, 20
        n_x_to, n_theta_to = 20, 18
        length = 30.
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=.2))
        if gui:
            case.plot()

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.13
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


if __name__ == '__main__':
    KratosUnittest.main()
