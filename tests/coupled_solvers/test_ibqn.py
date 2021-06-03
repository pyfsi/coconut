from coconut.tests.coupled_solvers import coupled_solver

import unittest


class TestCoupledSolverIBQN(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_ibqn.json'


if __name__ == '__main__':
    unittest.main()
