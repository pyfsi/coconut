from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
from coconut import tools
from coconut.data_structure import Model, Interface

import numpy as np
from os.path import join
import json


def create(parameters):
    return SolverWrapperTubeRingmodel(parameters)


class SolverWrapperTubeRingmodel(SolverWrapper):
    check_coupling_convergence_possible = True  # can solver check convergence after 1 iteration?

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

        # reading
        self.parameters = parameters
        self.settings = parameters['settings']
        self.working_directory = self.settings['working_directory']
        input_file = self.settings.get('input_file')
        if input_file is not None:
            case_file_name = join(self.working_directory, input_file)
            with open(case_file_name, 'r') as case_file:
                case_file_settings = json.load(case_file)
            case_file_settings.update(self.settings)
            with open(case_file_name, 'w') as case_file:
                json.dump(case_file_settings, case_file, indent=2)
            self.settings.update(case_file_settings)

        # settings
        l = self.settings['l']  # length
        d = self.settings['d']  # diameter
        self.rreference = d / 2  # reference radius of cross section
        self.rhof = self.settings['rhof']  # fluid density

        self.preference = self.settings.get('preference', 0)  # reference pressure

        e = self.settings['e']  # Young's modulus of structure
        h = self.settings['h']  # thickness of structure
        self.cmk2 = (e * h) / (self.rhof * d)  # wave speed squared

        self.m = self.settings['m']  # number of segments
        self.dz = l / self.m  # segment length
        axial_offset = self.settings.get('axial_offset', 0)  # start position along axis
        self.z = axial_offset + np.arange(self.dz / 2 - l / 2, l / 2, self.dz)  # data is stored in cell centers

        self.k = 0  # iteration
        self.n = 0  # time step

        self.residual_atol = self.settings.get('residual_atol')  # absolute residual tolerance for coupling convergence

        # initialization
        self.areference = np.pi * d ** 2 / 4  # reference area of cross section
        self.p = np.zeros(self.m) * self.preference  # kinematic pressure
        self.a = np.zeros(self.m) * self.areference  # area of cross section
        self.c02 = self.cmk2 - self.preference / 2  # wave speed squared with reference pressure

        self.disp = np.zeros((self.m, 3))  # displacement
        self.trac = np.zeros((self.m, 3))  # traction (always zero)

        # create ModelParts
        self.model = Model()
        self.input_model_part_name = self.settings['interface_input'][0]['model_part']
        self.output_model_part_name = self.settings['interface_output'][0]['model_part']
        self.model_part = self.model.create_model_part(self.input_model_part_name, np.zeros(self.m),
                                                       self.rreference * np.ones(self.m), self.z, np.arange(self.m))
        if self.input_model_part_name != self.output_model_part_name:
            self.model_part = self.model.create_model_part(self.output_model_part_name, np.zeros(self.m),
                                                           self.rreference * np.ones(self.m), self.z, np.arange(self.m))

        # interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_input.set_variable_data(self.input_model_part_name, 'pressure', self.p.reshape(-1, 1))
        self.interface_input.set_variable_data(self.input_model_part_name, 'traction', self.trac)
        self.interface_output = Interface(self.settings['interface_output'], self.model)
        self.interface_output.set_variable_data(self.output_model_part_name, 'displacement', self.disp)

        self.output_solution_step()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

        if self.check_coupling_convergence and self.residual_atol is None:
            raise ValueError(f'To check the coupling convergence with {self.__class__.__name__},'
                             f' the parameter "residual_atol" needs to be specified')

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.k = 0
        self.n += 1

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # input
        self.interface_input = interface_input.copy()
        self.p = interface_input.get_variable_data(self.input_model_part_name,
                                                   'pressure').flatten() / self.rhof  # kinematic pressure

        # coupling convergence
        if self.check_coupling_convergence:
            residual = np.linalg.norm(self.a - self.areference * (2 / (2 + (self.preference - self.p) / self.c02)) ** 2)
            if residual < self.residual_atol:
                self.coupling_convergence = True
                if self.print_coupling_convergence:
                    tools.print_info(f'{self.__class__.__name__} converged')

        # independent rings model
        for i in range(len(self.p)):
            if self.p[i] > 2 * self.c02 + self.preference:
                raise ValueError('Unphysical pressure')
        self.a = self.areference * (2 / (2 + (self.preference - self.p) / self.c02)) ** 2

        self.k += 1
        if self.debug:
            p = self.p * self.rhof
            file_name = join(self.working_directory, f'input_pressure_traction_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'x-traction':<22}"
                           f"\t{'y-traction':<22}\t{'z-traction':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{p[i]:<22}\t{self.trac[i, 0]:<22}'
                               f'\t{self.trac[i, 1]:<22}\t{self.trac[i, 2]:<22}\n')
            file_name = join(self.working_directory, f'output_area_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        # output
        self.disp[:, 1] = np.sqrt(self.a / np.pi) - self.rreference
        self.interface_output.set_variable_data(self.output_model_part_name, 'displacement', self.disp)
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        if self.debug:
            file_name = join(self.working_directory, f'area_ts{self.n}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

    def finalize(self):
        super().finalize()
