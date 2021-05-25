from coconut.tests.coupled_solvers import coupled_solver

import unittest


class TestCoupledSolverIQNI(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_iqni.json'


if __name__ == '__main__':
    unittest.main()
