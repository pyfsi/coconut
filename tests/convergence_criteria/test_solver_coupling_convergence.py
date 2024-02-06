from coconut.tools import create_instance
from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper

import unittest
from unittest.mock import Mock


class TestConvergenceCriterionSolverCouplingConvergence(unittest.TestCase):

    def test_initialize(self):
        # create convergence criterion
        parameters = {'type': 'convergence_criteria.solver_coupling_convergence',
                      'settings': {'solver_index': 0}}
        convergence_criterion = create_instance(parameters)

        with self.assertRaises(ValueError):
            convergence_criterion.initialize(Mock())

        convergence_criterion = create_instance(parameters)
        with self.assertRaises(ValueError):
            convergence_criterion.initialize([1, 2])

        convergence_criterion = create_instance(parameters)
        solver_wrappers = [Mock(spec=SolverWrapper), Mock(spec=SolverWrapper)]
        solver_wrappers[0].check_coupling_convergence_possible = False

        self.assertFalse(solver_wrappers[0].check_coupling_convergence_possible)
        with self.assertRaises(ValueError):
            convergence_criterion.initialize(solver_wrappers)

        convergence_criterion = create_instance(parameters)
        solver_wrappers[0].check_coupling_convergence_possible = True

        convergence_criterion.initialize(solver_wrappers)

        self.assertTrue(convergence_criterion.solver_wrapper.check_coupling_convergence)

    def test_is_satisfied(self):
        # test convergence criterion

        # create convergence criterion
        parameters = {'type': 'convergence_criteria.solver_coupling_convergence',
                      'settings': {'solver_index': 0}}
        convergence_criterion = create_instance(parameters)

        solver_wrappers = [Mock(spec=SolverWrapper), Mock(spec=SolverWrapper)]
        solver_wrappers[0].check_coupling_convergence_possible = True
        convergence_criterion.initialize(solver_wrappers)

        solver_wrappers[0].coupling_convergence = False
        self.assertFalse(convergence_criterion.is_satisfied())

        solver_wrappers[0].coupling_convergence = True
        self.assertTrue(convergence_criterion.is_satisfied())


if __name__ == '__main__':
    unittest.main()
