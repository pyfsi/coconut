from coconut.tools import create_instance, get_solver_env, solver_available, print_box
import coconut.coupling_components.solver_wrappers.openfoam.openfoam_io as of_io
from coconut.data_structure import Model
from coconut.data_structure import Interface

import unittest
import numpy as np
import os
from os.path import join
import multiprocessing
import shutil
import re
import json
from subprocess import check_call, DEVNULL


class TestSolverWrapperOpenFOAMExtend(unittest.TestCase):
    version = 41  # OpenFOAM version without dot, e.g. 41 , set in sub-class

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to openfoam directory
        cls.file_name = join(dir_name, f'test_v{cls.version}_updated_solver/wire/parameters.json')
        cls.working_dir = join(dir_name, f'test_v{cls.version}/wire/CSM')

        # setup
        shutil.rmtree(os.path.join(dir_name, cls.working_dir), ignore_errors=True)
        shutil.copytree(os.path.join(dir_name, f'test_v{cls.version}/wire/setup'), cls.working_dir)

    def setUp(self):
        with open(self.file_name, 'r') as parameter_file:
            self.parameters = json.load(parameter_file)
        self.parameters['settings']['working_directory'] = os.path.relpath(self.working_dir)  # set working directory

        self.mp_name_in = self.parameters['settings']['interface_input'][0]['model_part']
        self.mp_name_out = self.parameters['settings']['interface_output'][0]['model_part']

        self.folder_path = os.path.join(os.getcwd(), self.parameters['settings']['working_directory'])
        self.delta_t = self.parameters['settings']['delta_t']
        self.t_prec = self.parameters['settings']['time_precision']
        self.max_cores = min(4, multiprocessing.cpu_count())  # number of cores for parallel calculation

        solver_name = self.parameters['type'].replace('solver_wrappers.', '')
        self.env = get_solver_env(solver_name, self.folder_path)
        # self.clean_case()
        self.set_up_case()


    # def tearDown(self):
    #     self.clean_case()
    #
    # @classmethod
    # def tearDownClass(cls):
    #     shutil.rmtree(cls.working_dir)

    def clean_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allclean'), shell=True, env=self.env)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'prepareCase'), shell=True, env=self.env)

    def set_cores(self, cores):
        if cores == 1:
            self.parameters['settings']['parallel'] = False
        else:
            self.parameters['settings']['parallel'] = True
            decompose_file_name = os.path.join(self.folder_path, 'system', 'decomposeParDict')
            with open(decompose_file_name, 'r') as f:
                dct = f.read()
            new_dict = re.sub(r'numberOfSubdomains[\s\n]+' + of_io.int_pattern, f'numberOfSubdomains   {cores}', dct)

            with open(decompose_file_name, 'w') as f:
                f.write(new_dict)

    def get_p(self, x):
        return np.sin(x) + 9e8
        # return x * 0

    def get_p0(self, x):
        return np.sin(x) + 1e5
        # return x * 0

    def get_shear(self, x):
        shear = np.zeros((x.shape[0], 3))
        return shear

    def test_displacement(self):
        # adapt parameters, create solver
        self.set_cores(1)
        solver = create_instance(self.parameters)
        solver.initialize()
        timesteps = 2

        # set number of time steps.
        for i in range(timesteps):

            solver.initialize_solution_step()
            interface_input = solver.get_interface_input()

            #set pressure and traction
            model_part=interface_input.get_model_part(self.mp_name_in)
            x0 = model_part.x0
            p = self.get_p(x0)
            p0 = self.get_p0(x0)
            t = self.get_shear(x0)
            pr = np.column_stack(p)
            pr0 = np.column_stack(p0)
            tr = np.column_stack(t)
            pr_list = [pr,np.zeros_like(pr),pr]
            tr_list = [np.zeros_like(tr)]

            #run solver for three pressure and traction(first one = last one)
            displacement = []

            for pr in pr_list:
                pressure = pr.reshape((-1, 1))
                interface_input.set_variable_data(self.mp_name_in, "pressure", pressure)
                interface_output = solver.solve_solution_step(interface_input)
                displacement.append(interface_output.get_variable_data(self.mp_name_out,'displacement'))
            solver.finalize_solution_step()

            print(displacement[0])
            print(displacement[1])
            print(displacement[2])

            # check if same position give same displacement
            np.testing.assert_allclose(displacement[0], displacement[2], atol=1e-9, rtol=0)

        solver.finalize()

if __name__ == "__main__":
    unittest.main()




