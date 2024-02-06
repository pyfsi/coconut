from coconut.tools import create_instance, cd
from coconut.tests.coupled_solvers import coupled_solver

import unittest
import numpy as np


class TestCoupledSolverExplicit(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_explicit.json'

    def test_coupled_solver(self):
        with cd(self.working_dir):
            coupled_solver = create_instance(self.parameters)
            coupled_solver.initialize()

            coupled_solver.initialize_solution_step()
            coupled_solver.solve_solution_step()

            sol_x = [0.00000000e+00, 3.91249330e-07, 0.00000000e+00, 0.00000000e+00, 3.80784228e-07, 0.00000000e+00,
                     0.00000000e+00, 3.70319125e-07, 0.00000000e+00, 0.00000000e+00, 3.59854023e-07, 0.00000000e+00,
                     0.00000000e+00, 3.49388922e-07, 0.00000000e+00, 0.00000000e+00, 3.38923821e-07, 0.00000000e+00,
                     0.00000000e+00, 3.28458720e-07, 0.00000000e+00, 0.00000000e+00, 3.17993620e-07, 0.00000000e+00,
                     0.00000000e+00, 3.07528521e-07, 0.00000000e+00, 0.00000000e+00, 2.97063421e-07, 0.00000000e+00]

            # TODO: The reference solution has to be modified. Future work: mock solver

            np.testing.assert_allclose(coupled_solver.x.get_interface_data(), sol_x, rtol=1e-5)

            coupled_solver.finalize_solution_step()
            coupled_solver.output_solution_step()

            coupled_solver.finalize()


if __name__ == '__main__':
    unittest.main()
