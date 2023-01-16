from coconut.tools import create_instance, cd
from coconut.tests.coupled_solvers import coupled_solver

import unittest
import numpy as np


class TestCoupledSolverOneway(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_one_way.json'

    def test_coupled_solver(self):
        with cd(self.working_dir):
            coupled_solver = create_instance(self.parameters)
            coupled_solver.initialize()

            coupled_solver.initialize_solution_step()
            coupled_solver.solve_solution_step()

            sol_y = [46.9498829, 45.69407251, 44.43826212, 43.18245173, 41.92664134,
                     40.67083095, 39.41502056, 38.15921017, 36.90339978, 35.64758939,
                     0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0.]

            # TODO: The reference solution has to be modified. Future work: mock solver

            np.testing.assert_allclose(coupled_solver.y.get_interface_data(), sol_y, rtol=1e-5)

            coupled_solver.finalize_solution_step()
            coupled_solver.output_solution_step()

            coupled_solver.finalize()


if __name__ == '__main__':
    unittest.main()
