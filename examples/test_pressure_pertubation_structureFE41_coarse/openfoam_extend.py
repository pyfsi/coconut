from coconut.tools import create_instance, get_solver_env, solver_available, print_box
import coconut.coupling_components.solver_wrappers.openfoam.openfoam_io as of_io
import matplotlib.pyplot as plt
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
        cls.file_name = join(dir_name, f'test_v{cls.version}/wire/parameters.json')
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
        self.iterations = self.parameters['settings']['iterations']
        self.boundary_names = self.parameters['settings']['boundary_names']
        self.max_cores = min(4, multiprocessing.cpu_count())  # number of cores for parallel calculation

        solver_name = self.parameters['type'].replace('solver_wrappers.', '')
        self.env = get_solver_env(solver_name, self.folder_path)
        # self.clean_case()
        self.set_up_case()

    def clean_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'Allclean'), shell=True, env=self.env)

    def set_up_case(self):
        check_call('sh ' + os.path.join(self.folder_path, 'prepareCase'), shell=True, env=self.env)

    def read_face_centres(self, boundary_name, nfaces):
        filename_ccx = os.path.join("/lusers/temp/mathieu/PycharmProjects/coconut/examples/test_pressure_pertubation_structureFE41_coarse/test_v41/wire/CSM",'0.00010','ccx')
        filename_ccy = os.path.join("/lusers/temp/mathieu/PycharmProjects/coconut/examples/test_pressure_pertubation_structureFE41_coarse/test_v41/wire/CSM",'0.00010','ccy')
        filename_ccz = os.path.join("/lusers/temp/mathieu/PycharmProjects/coconut/examples/test_pressure_pertubation_structureFE41_coarse/test_v41/wire/CSM",'0.00010','ccz')

        x = of_io.get_boundary_field(file_name=filename_ccx, boundary_name=boundary_name, size=nfaces,
                                     is_scalar=True)
        y = of_io.get_boundary_field(file_name=filename_ccy, boundary_name=boundary_name, size=nfaces,
                                     is_scalar=True)
        z = of_io.get_boundary_field(file_name=filename_ccz, boundary_name=boundary_name, size=nfaces,
                                     is_scalar=True)
        return x, y, z

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

    def test_displacement(self):
        # adapt parameters, create solver
        self.set_cores(1)
        solver = create_instance(self.parameters)
        solver.initialize()
        timesteps = 1

        # set number of time steps.
        for i in range(timesteps):

            solver.initialize_solution_step()
            interface_input = solver.get_interface_input()

            #set pressure and traction
            model_part=interface_input.get_model_part(self.mp_name_in)
            x0 = model_part.x0
            p = np.array([950731.762085373, 950731.762085373, 734882379.7858711, 1264865992.9069302, 1468537141.5581686, 114910630.44785258])
            p = p*1.5
            print(p)
            pr = np.column_stack(p)
            pp = np.zeros(pr.size)

            pr_list = []
            displacement = []
            displacement_field = []
            cell_centres = []

            for i in range(self.iterations):
                if i == 0:
                    pp = pr
                elif i == self.iterations-1:
                    pp = pr
                else:
                    pp = pr + 0.01 * pr / (2 ** (i - 1))

                pr_list.append(pp)

            for pr in pr_list:
                pressure = pr.reshape((-1, 1))
                interface_input.set_variable_data(self.mp_name_in, "pressure", pressure)
                interface_output = solver.solve_solution_step(interface_input)

                mp = interface_output.get_model_part(self.mp_name_out)
                nfaces = mp.size
                filename_displacement = os.path.join("/lusers/temp/mathieu/PycharmProjects/coconut/examples/test_pressure_pertubation_structureFE41_coarse/test_v41/wire/CSM",'0.00010','U')
                disp_field = of_io.get_boundary_field(file_name=filename_displacement, boundary_name="wireTopContact",
                                                      size=nfaces, is_scalar=False)
                x, y, z = self.read_face_centres("wireTopContact", nfaces)
                rad_disp = disp_field[:,1]
                cell_centres.append(x)
                displacement_field.append(rad_disp)

                displacement.append(interface_output.get_variable_data(self.mp_name_out,'displacement'))
            solver.finalize_solution_step()

            # fig, axs =plt.subplots(3)
            # fig.suptitle('Vertically stacked subplots')
            # for j in range(self.iterations):
            #     diff = displacement_field[j] - displacement_field[0]
            #     axs[0].plot(cell_centres[j],diff, label = f'iteration_{j}')
            #     axs[0].legend(loc = 2)
            #     axs[1].plot(cell_centres[j],displacement_field[j], label = f'displacement_field_{j}')
            #     axs[1].legend(loc = 2)
            #     k = pr_list[j]
            #     axs[2].plot(x0,k[0,:], label = f'pressure_distribution_{j}')
            #     axs[2].legend(loc = 2)
            #
            plt.legend()
            plt.show()

            # print("displacement")
            difference = displacement_field[4]-displacement_field[0]
            plt.plot(difference)
            plt.show()
            print("cell_centre_diff")
            n =  cell_centres[4]-cell_centres[0]
            for u in range(n.size):
                print('%.18f' % n[u])
            v = pr_list[0]
            w = pr_list[4]
            print("pressure_diff")
            z =v[0,:]-w[0,:]
            for i in range(z.size):
                print('%.18f' % z[i])

            # check if same position give same displacement
            np.testing.assert_allclose(displacement[0], displacement[0], atol=1e-14, rtol=0)

        solver.finalize()

if __name__ == "__main__":
    unittest.main()




