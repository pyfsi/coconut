from unittest import TestSuite, TestLoader, TestCase
import unittest
import sys


def key_in_test(key, test):
    test_suite = TestSuite()
    if isinstance(test, TestSuite):
        if test.countTestCases() == 0:
            return []
        for test_case in test:
            test_suite.addTests(key_in_test(key, test_case))
        return test_suite
    elif isinstance(test, TestCase):
        if key.lower() in str(test).lower():
            test_suite.addTest(test)
    else:
        raise TypeError('Variable "test" is not an instance of TestSuite or TestCase')
    return test_suite


def generate_test_suite(keys, all_tests):
    if not keys:  # empty list
        return generate_test_suite(['-fast'], all_tests)
    test_suite = TestSuite()
    for key in keys:
        if key == '-all':
            test_suite.addTests(all_tests)
            return test_suite
        elif key == '-fast':
            test_suite.addTests(generate_test_suite(['tests.convergence_criteria', 'tests.coupled_solvers',
                                                     'tests.data_structure', 'tests.mappers', 'tests.predictors',
                                                     'tests.solver_wrappers.python', 'tests.solver_wrappers.combined',
                                                     'tests.post_processing'],
                                                    all_tests))
        else:
            for test in all_tests:
                test_suite.addTests(key_in_test(key, test))
    return test_suite


if __name__ == '__main__':
    keys = sys.argv[1:]
    loader = TestLoader()
    all_tests = loader.discover('coconut.tests')
    runner = unittest.TextTestRunner(verbosity=2, buffer=True)
    runner.run(generate_test_suite(keys, all_tests))
