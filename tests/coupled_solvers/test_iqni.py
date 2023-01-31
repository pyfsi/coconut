from coconut.tests.coupled_solvers import coupled_solver

import unittest


class TestCoupledSolverIQNI(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_iqni.json'

    def set_new_values(self):
        self.omega_new = 0.1
        self.settings['omega'] = self.omega_new

    def check_new_values(self, coupled_solver):
        self.assertEqual(coupled_solver.omega, self.omega_new)


if __name__ == '__main__':
    unittest.main()
