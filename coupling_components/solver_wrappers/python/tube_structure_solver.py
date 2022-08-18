from coconut.coupling_components.component import Component
from coconut import tools
from coconut.data_structure import Model, Interface
import coconut.coupling_components.solver_wrappers.python.banded as bnd

import numpy as np
import os
from os.path import join
from scipy.linalg import solve_banded
import json
import pickle


def create(parameters):
    return SolverWrapperTubeStructure(parameters)


class SolverWrapperTubeStructure(Component):
    al = 2  # Number of terms below diagonal in matrix
    au = 2  # Number of terms above diagonal in matrix

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__()

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
        self.unsteady = self.settings.get('unsteady', True)
        self.solver = self.settings.get('solver', 'solve_banded')  # method for solving the linear system
        if self.solver not in ('solve_banded', 'direct'):
            raise ValueError('The value of key "solver" must be "direct" or "solve_banded"')
        self.timestep_start = self.settings.get('timestep_start', 0)
        self.save_restart = self.settings.get('save_restart', 0)  # eg to restart

        l = self.settings['l']  # length
        d = self.settings['d']  # diameter
        self.rreference = d / 2  # reference radius of cross section
        self.rhof = self.settings['rhof']  # fluid density

        self.preference = self.settings.get('preference', 0)  # reference pressure

        e = self.settings['e']  # Young's modulus of structure
        nu = self.settings['nu']  # Poisson's ratio
        self.h = self.settings['h']  # thickness of structure
        self.rhos = self.settings['rhos']  # structure density
        self.cmk2 = (e * self.h) / (self.rhof * d)  # wave speed squared
        self.b1 = (self.h * e) / (1 - nu ** 2) * (self.h ** 2) / 12
        self.b2 = self.b1 * (2 * nu) / self.rreference ** 2
        self.b3 = (self.h * e) / (1 - nu ** 2) * 1 / self.rreference ** 2

        self.m = self.settings['m']  # number of segments
        self.dz = l / self.m  # segment length
        axial_offset = self.settings.get('axial_offset', 0)  # start position along axis
        self.z = axial_offset + np.arange(self.dz / 2 - l / 2, l / 2, self.dz)  # data is stored in cell centers

        self.k = 0  # iteration
        self.n = self.timestep_start  # time step
        if self.unsteady:
            self.dt = self.settings['delta_t']  # time step size
            self.time_discretization = self.settings.get('time_discretization', 'backward euler').lower()
            if self.time_discretization == 'newmark':
                self.gamma = self.settings['gamma']  # Newmark parameter: gamma >= 1/2
                self.beta = self.settings['beta']  # Newmark parameter: beta >= 1/4 * (1/2 + gamma) ^ 2
                self.nm = True
                if not self.gamma >= 0.5 or not self.beta >= 0.25 * (0.5 + self.gamma) ** 2:
                    raise Exception('Inadequate Newmark parameteres')
            elif self.time_discretization == 'backward euler':
                self.gamma = 1
                self.beta = 1
                self.nm = False  # used to set rnddot to zero
            else:
                raise ValueError('Time discretization should be \"Newmark\" or \"backward Euler\"')
        else:
            # not used
            self.dt = 1  # Time step size default
            self.beta = 1
            self.nm = False

        # initialization
        self.areference = np.pi * self.rreference ** 2  # reference area of cross section
        self.p = np.ones(self.m) * self.preference  # pressure
        self.a = np.ones(self.m) * self.areference  # area of cross section
        if self.timestep_start == 0:  # no restart
            self.r = np.ones(self.m + 4) * self.rreference  # radius of cross section
            if self.unsteady:
                self.rdot = np.zeros(self.m)  # first derivative of the radius with respect to time in current timestep
                self.rddot = np.zeros(self.m)  # second derivative of radius with respect to time in current timestep
        else:  # restart
            file_name = join(self.working_directory, f'case_timestep{self.timestep_start}.pickle')
            with open(file_name, 'rb') as file:
                data = pickle.load(file)
            self.r = data['r']  # radius of cross section
            self.rdot = data['rdot']  # first derivative of the radius with respect to time in current timestep
            self.rddot = data['rddot']  # second derivative of the radius with respect to time in current timestep
        self.rn = np.array(self.r)  # previous radius of cross section
        self.rndot = np.zeros(self.m)  # first derivative of the radius with respect to time in previous timestep
        self.rnddot = np.zeros(self.m)  # second derivative of the radius with respect to time in previous timestep

        self.disp = np.zeros((self.m, 3))  # displacement
        self.trac = np.zeros((self.m, 3))  # traction (always zero)

        self.conditioning = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.unsteady
                             + 6.0 * self.b1 / self.dz ** 4
                             + 2.0 * self.b2 / self.dz ** 2 + self.b3)  # factor for conditioning Jacobian

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

        # calculate inverse Jacobian if direct solver
        if self.solver == 'direct':
            j = bnd.to_dense(self.get_jacobian())
            self.ji = np.linalg.inv(j)

        # time
        self.init_time = self.init_time
        self.run_time = 0.0

        # debug
        self.debug = self.settings.get('debug', False)  # save input and output of each iteration of every time step
        self.output_solution_step()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.k = 0
        self.n += 1
        if self.unsteady:
            self.rn = np.array(self.r)
            self.rndot = self.rdot
            self.rnddot = self.rddot

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # input
        self.interface_input = interface_input.copy()
        self.p = interface_input.get_variable_data(self.input_model_part_name, 'pressure').flatten()

        # solve system
        if self.solver == 'solve_banded':
            j = self.get_jacobian()
            b = self.get_b()
            self.r = solve_banded((self.al, self.au), j, b)
        elif self.solver == 'direct':
            b = self.get_b()
            self.r = self.ji @ b

        self.k += 1
        if self.debug:
            file_name = join(self.working_directory, f'input_pressure_traction_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'x-traction':<22}"
                           f"\t{'y-traction':<22}\t{'z-traction':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.p[i]:<22}\t{self.trac[i, 0]:<22}'
                               f'\t{self.trac[i, 1]:<22}\t{self.trac[i, 2]:<22}\n')
            file_name = join(self.working_directory, f'output_area_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        # output does not contain boundary conditions
        self.a = self.r[2:self.m + 2] ** 2 * np.pi
        self.disp[:, 1] = self.r[2:self.m + 2] - self.rreference
        self.interface_output.set_variable_data(self.output_model_part_name, 'displacement', self.disp)
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

        if self.unsteady:
            self.rddot = ((self.r[2:self.m + 2] - self.rn[2:self.m + 2]) / (self.beta * self.dt ** 2)
                          - self.rndot / (self.beta * self.dt) - self.rnddot * (1 / (2 * self.beta) - 1) * self.nm)
            self.rdot = self.rndot + self.dt * (1 - self.gamma) * self.rnddot + self.dt * self.gamma * self.rddot

    def finalize(self):
        super().finalize()

    def output_solution_step(self):
        if self.n > 0 and self.save_restart != 0 and self.n % self.save_restart == 0:
            file_name = join(self.working_directory, f'case_timestep{self.n}.pickle')
            with open(file_name, 'wb') as file:
                pickle.dump({'r': self.r, 'rdot': self.rdot, 'rddot': self.rddot}, file)
            if self.save_restart < 0 and self.n + self.save_restart > self.timestep_start:
                try:
                    os.remove(join(self.working_directory, f'case_timestep{self.n + self.save_restart}.pickle'))
                except OSError:
                    pass
        if self.debug:
            file_name = join(self.working_directory, f'area_ts{self.n}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

    def get_interface_input(self):
        return self.interface_input.copy()

    def get_interface_output(self):
        return self.interface_output.copy()

    def get_residual(self):
        f = np.zeros(self.m + 4)
        f[0] = (self.r[0] - self.rreference)
        f[1] = (self.r[1] - self.rreference)
        f[2:self.m + 2] = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.r[2: self.m + 2] * self.unsteady
                           + self.b1 / self.dz ** 4 * (self.r[4:self.m + 4] - 4.0 * self.r[3:self.m + 3]
                                                       + 6.0 * self.r[2:self.m + 2] - 4.0 * self.r[1:self.m + 1]
                                                       + self.r[0:self.m])
                           - self.b2 / self.dz ** 2 * (self.r[3:self.m + 3] - 2.0 * self.r[2:self.m + 2]
                                                       + self.r[1:self.m + 1])
                           + self.b3 * (self.r[2:self.m + 2] - self.rreference)
                           - (self.p - self.preference)
                           - self.rhos * self.h * (self.rn[2:self.m + 2] / (self.beta * self.dt ** 2)
                                                   + self.rndot / (self.beta * self.dt)
                                                   + self.rnddot * (1.0 / (2.0 * self.beta) - 1.0) * self.nm)
                           * self.unsteady) / self.conditioning
        f[self.m + 2] = (self.r[self.m + 2] - self.rreference)
        f[self.m + 3] = (self.r[self.m + 3] - self.rreference)
        return f

    def get_b(self):
        f = np.zeros(self.m + 4)
        f[0] = self.rreference
        f[1] = self.rreference
        f[2:self.m + 2] = (self.b3 * self.rreference
                           + (self.p - self.preference)
                           + self.rhos * self.h * (self.rn[2:self.m + 2] / (self.beta * self.dt ** 2)
                                                   + self.rndot / (self.beta * self.dt)
                                                   + self.rnddot * (1.0 / (2.0 * self.beta) - 1.0) * self.nm)
                           * self.unsteady) / self.conditioning
        f[self.m + 2] = self.rreference
        f[self.m + 3] = self.rreference
        return f

    def get_jacobian(self):
        j = np.zeros((self.al + self.au + 1, self.m + 4))
        j[self.au + 0 - 0, 0] = 1.0  # [0, 0]
        j[self.au + 1 - 1, 1] = 1.0  # [1, 1]
        j[self.au + 2, 0:self.m] = self.b1 / self.dz ** 4 / self.conditioning  # [i, (i - 2)]
        j[self.au + 1, 1:self.m + 1] = (-4.0 * self.b1 / self.dz ** 4
                                        - self.b2 / self.dz ** 2) / self.conditioning  # [i, (i - 1)]
        j[self.au + 0, 2:self.m + 2] = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.unsteady
                                        + 6.0 * self.b1 / self.dz ** 4 + 2.0 * self.b2 / self.dz ** 2
                                        + self.b3) / self.conditioning  # [i, i]
        j[self.au - 1, 3:self.m + 3] = (-4.0 * self.b1 / self.dz ** 4
                                        - self.b2 / self.dz ** 2) / self.conditioning  # [i, (i + 1)]
        j[self.au - 2, 4:self.m + 4] = self.b1 / self.dz ** 4 / self.conditioning  # [i, (i + 2)]
        j[self.au + (self.m + 2) - (self.m + 2), self.m + 2] = 1.0  # [m + 2, m + 2]
        j[self.au + (self.m + 3) - (self.m + 3), self.m + 3] = 1.0  # [m + 3, m + 3]
        return j

    def get_surrogate_jacobian(self):
        # df/dr
        js = self.get_jacobian()
        js = bnd.remove_boundaries(js, 2) * self.conditioning
        js = bnd.to_dense(js)

        # dr/df
        jis = np.linalg.inv(js)

        # df/dp
        jp = -1

        # dr/dp
        jsurrogate = jis * -jp
        return jsurrogate
