from coconut import data_structure
from coconut.data_structure.interface import Interface
from coconut.tools import create_instance, cd
from coconut.tests.coupled_solvers import coupled_solver

import unittest
import numpy as np
from coconut.tools import print_info

class TestCoupledSolverOneway(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_oneway.json'

    def test_coupled_solver(self):
        with cd(self.working_dir):
            coupled_solver = create_instance(self.parameters)
            coupled_solver.initialize()

            coupled_solver.initialize_solution_step()
            coupled_solver.solve_solution_step()

            sol_x = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

            # TODO: The reference solution has to be modified. Future work: mock solver

            np.testing.assert_allclose(coupled_solver.x.get_interface_data(), sol_x, rtol=1e-5)

            coupled_solver.finalize_solution_step()
            coupled_solver.output_solution_step()

            coupled_solver.finalize()

if __name__ == '__main__':
    unittest.main()
