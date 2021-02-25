from coconut.tools import create_instance

import unittest
import os
import json
import numpy as np


class TestCoupledSolverRelaxation(unittest.TestCase):

    def test_coupled_solver_relaxation(self):
        # read settings
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_relaxation.json')
        with open(parameter_file_name, 'r') as parameter_file:
            settings = json.load(parameter_file)

        coupled_solver = create_instance(settings)
        coupled_solver.initialize()

        coupled_solver.initialize_solution_step()
        coupled_solver.solve_solution_step()
        sol_x = [0.00000e+00, 3.09851e-07, 0.00000e+00, 0.00000e+00,
                 3.00094e-07, 0.00000e+00, 0.00000e+00, 2.90572e-07,
                 0.00000e+00, 0.00000e+00, 2.81238e-07, 0.00000e+00,
                 0.00000e+00, 2.72094e-07, 0.00000e+00, 0.00000e+00,
                 2.63131e-07, 0.00000e+00, 0.00000e+00, 2.54343e-07,
                 0.00000e+00, 0.00000e+00, 2.45726e-07, 0.00000e+00,
                 0.00000e+00, 2.37273e-07, 0.00000e+00, 0.00000e+00,
                 2.28979e-07, 0.00000e+00]
        # TODO: The reference solution has to be modified. Future work: mock solver
        np.testing.assert_allclose(coupled_solver.x.get_interface_data(), sol_x, rtol=1e-1)

        coupled_solver.finalize_solution_step()
        coupled_solver.output_solution_step()

        coupled_solver.finalize()


if __name__ == '__main__':
    unittest.main()
