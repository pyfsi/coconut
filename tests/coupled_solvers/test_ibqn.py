from coconut.tests.coupled_solvers import coupled_solver

import unittest


class TestCoupledSolverIBQN(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_ibqn.json'

    def set_new_values(self):
        self.omega_new = 0.1
        self.atol_new = 1e-10
        self.rtol_new = 1e-15
        self.settings['omega'] = self.omega_new
        self.settings['absolute_tolerance_gmres'] = self.atol_new
        self.settings['relative_tolerance_gmres'] = self.rtol_new

    def check_new_values(self, coupled_solver):
        self.assertEqual(coupled_solver.omega, self.omega_new)
        self.assertEqual(coupled_solver.atol, self.atol_new)
        self.assertEqual(coupled_solver.rtol, self.rtol_new)


if __name__ == '__main__':
    unittest.main()
