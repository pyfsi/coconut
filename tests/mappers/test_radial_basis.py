from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.tests.mappers.test_interpolator import Case1D, Case2D, Case3DSphere, Case3DCylinder, Case3DSinc
from coconut.coupling_components.tools import quicktimer

import os


class TestMapperRadialBasis(KratosUnittest.TestCase):

    def test_mapper_radial_basis(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_radial_basis.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_mapper = parameters['mapper']

        gui = 0  # *** gui gives problems when running all tests?

        # 1D case: square-root grid + linear function
        """
        n_from = 14, n_to = 5 
            => max error = 5.8e-5
        """
        n_from, n_to = 14, 5
        par_mapper['settings'].SetArray('directions', ['Z'])

        case = Case1D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=1e-4))
        if gui:
            case.plot()

        # 2D case: circle + linear function
        """
        n_from = 33, n_to = 22 
            => max error = 0.003
        """
        n_from, n_to = 33, 22
        par_mapper['settings'].SetArray('directions', ['X', 'Y'])

        case = Case2D(n_from, n_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=0.005))
        if gui:
            case.plot()

        # 3D case: sphere + sine function
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
            => max error = 3.2e-4
        """
        n_theta_from, n_phi_from = 50, 30
        n_theta_to, n_phi_to = 22, 11
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSphere(n_theta_from, n_phi_from, n_theta_to, n_phi_to)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=5e-4))
        if gui:
            case.plot()

        # 3D case: cylinder + sine function
        """
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
            => max error = 0.00047
        """
        n_x_from, n_theta_from = 50, 15
        n_x_to, n_theta_to = 45, 13
        length = 1000.
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
        par_mapper['settings'].AddEmptyValue('scaling')
        par_mapper['settings'].SetArray('scaling', [.01, 1, 1])  # bad result without scaling

        case = Case3DCylinder(n_x_from, n_theta_from, n_x_to, n_theta_to, length)
        case.map(par_mapper)
        self.assertTrue(case.check(tolerance=5e-4))
        if gui:
            case.plot()
        par_mapper['settings'].RemoveValue('scaling')

        # 3D case: sinc + linear vector function
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
            => max error = 0.05
            
        n_x_from, n_y_from = 50, 50
        n_x_to, n_y_to = 60, 60
            => max error = 0.0024
        """
        n_x_from, n_y_from = 20, 20
        n_x_to, n_y_to = 13, 13
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        case = Case3DSinc(n_x_from, n_y_from, n_x_to, n_y_to)
        case.map(par_mapper)
        for tmp in case.check(tolerance=0.1):
            self.assertTrue(tmp)
        if gui:
            case.plot()


if __name__ == '__main__':
    KratosUnittest.main()
