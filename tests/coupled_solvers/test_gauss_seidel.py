from coconut.tests.coupled_solvers import coupled_solver

import unittest


class TestCoupledSolverGaussSeidel(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_gauss_seidel.json'


if __name__ == '__main__':
    unittest.main()
