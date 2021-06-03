from coconut.coupling_components.component import Component
from coconut import tools
from coconut.data_structure import Model, Interface

import numpy as np
import os.path as path
import json


def create(parameters):
    return SolverWrapperTubeRingmodel(parameters)


class SolverWrapperTubeRingmodel(Component):
    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

        # reading
        self.parameters = parameters
        self.settings = parameters['settings']
        self.working_directory = self.settings['working_directory']
        input_file = self.settings['input_file']
        case_file_name = path.join(self.working_directory, input_file)
        with open(case_file_name, 'r') as case_file:
            self.settings.update(json.load(case_file))

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

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

        # debug
        self.debug = False  # set on True to save input and output of each iteration of every time step
        self.output_solution_step()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

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

        # independent rings model
        for i in range(len(self.p)):
            if self.p[i] > 2 * self.c02 + self.preference:
                raise ValueError('Unphysical pressure')
        self.a = self.areference * (2 / (2 + (self.preference - self.p) / self.c02)) ** 2

        self.k += 1
        if self.debug:
            p = self.p * self.rhof
            file_name = self.working_directory + f'/Input_Pressure_Traction_TS{self.n}_IT{self.k}.txt'
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'x-traction':<22}"
                           f"\t{'y-traction':<22}\t{'z-traction':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{p[i]:<22}\t{self.trac[i, 0]:<22}'
                               f'\t{self.trac[i, 1]:<22}\t{self.trac[i, 2]:<22}\n')
            file_name = self.working_directory + f'/Output_Area_TS{self.n}_IT{self.k}.txt'
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        # output
        self.disp[:, 1] = np.sqrt(self.a / np.pi) - self.rreference
        self.interface_output.set_variable_data(self.output_model_part_name, 'displacement', self.disp)
        return self.interface_output  # TODO: make copy?

    def finalize_solution_step(self):
        super().finalize_solution_step()

    def finalize(self):
        super().finalize()

    def output_solution_step(self):
        if self.debug:
            file_name = self.working_directory + f'/Area_TS{self.n}.txt'
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

    def get_interface_input(self):  # TODO: need to have latest data?
        return self.interface_input

    def get_interface_output(self):  # TODO: need to have latest data?
        return self.interface_output
