from coconut.tests.coupled_solvers import coupled_solver
from coconut.tools import remove_recursively

import unittest
import os
import shutil


class TestCoupledSolverIQNISM(coupled_solver.TestCoupledSolver):
    parameter_file_name = 'test_iqnism.json'

    def setUp(self):
        super().setUp()

        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to coupled_solvers directory

        # set working directories
        self.working_dir = os.path.join(dir_name, 'coupled_solver_tmp')
        working_dir_sur_cfd = os.path.join(self.working_dir, 'CFD_surrogate')
        working_dir_sur_csm = os.path.join(self.working_dir, 'CSM_surrogate')
        settings = self.parameters['settings']
        settings['surrogate']['coupled_solver']['solver_wrappers'][0]['settings']['working_directory'] \
            = os.path.relpath(working_dir_sur_cfd, start=self.working_dir)
        settings['surrogate']['coupled_solver']['solver_wrappers'][1]['settings']['working_directory'] \
            = os.path.relpath(working_dir_sur_csm, start=self.working_dir)

        # setup
        os.mkdir(working_dir_sur_cfd)
        os.mkdir(working_dir_sur_csm)
        shutil.copy(os.path.join(dir_name, 'setup_tube_flow/solver_parameters.json'), working_dir_sur_cfd)
        shutil.copy(os.path.join(dir_name, 'setup_tube_ringmodel/solver_parameters.json'), working_dir_sur_csm)

    def remove_keys(self, keys):
        # remove keys to avoid warnings
        for key in keys:
            remove_recursively(key, self.parameters['solver_wrappers'])
            remove_recursively(key, self.parameters['settings']['model'])
            remove_recursively(key, self.parameters['settings']['surrogate'])


if __name__ == '__main__':
    unittest.main()
