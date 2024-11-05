from coconut.coupling_components.solver_wrappers.solver_wrapper import SolverWrapper
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
    return SolverWrapperSaturatedSolid(parameters)


class SolverWrapperSaturatedSolid(SolverWrapper):
    check_coupling_convergence_possible = False  # can solver check convergence after 1 iteration?

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

        # store settings
        self.dt = self.settings['timestep_size']
        self.timestep = self.settings['timestep_start']
        self.interface_settings = self.settings['interface']
        self.material_properties = self.settings['material_properties']

        self.x0 = self.interface_settings['x0']
        self.y0 = self.interface_settings['y0']
        self.x1 = self.interface_settings['x1']
        self.x1 = self.interface_settings['y1']
        self.n = self.interface_settings['faces']

        self.rho = self.material_properties['rho']
        self.L = self.material_proprties['latent']

        # initialization
        if self.timestep_start == 0:  # no restart
            x = np.linspace(self.x0, self.x1, self.n + 1)
            y = np.linspace(self.y0, self.y1, self.n + 1)
            ini_coord = []
            for i in range(self.n):
                ini_coord.append([(x[i] + x[i+1]) / 2, (y[i] + y[i+1]) / 2])
            self.ini_coord = np.array(ini_coord)

            self.prev_disp = np.zeros((self.n, 2))  # previous time step displacement
            self.dx = np.zeros((self.n, 2))  # latest time step displacement
            self.heat_flux = np.zeros(self.n) # heat flux [W/m^2]

        else:  # restart
            file_name = join(self.working_directory, f'case_timestep{self.timestep}.pickle')
            with open(file_name, 'rb') as file:
                data = pickle.load(file)
            self.ini_coord = data['ini_coord']
            self.prev_disp = data['prev_disp']
            self.dx = data['dx']
            self.heat_flux = data['heat_flux']

        # create ModelParts
        self.model = Model()
        self.input_model_part_name = self.settings['interface_input'][0]['model_part']
        self.output_model_part_name = self.settings['interface_output'][0]['model_part']
        self.model_part = self.model.create_model_part(self.input_model_part_name, self.ini_coord[:,0].flatten(),
                                                       self.ini_coord[:,1].flatten(), np.zeros(self.n), np.arange(self.n))
        if self.input_model_part_name != self.output_model_part_name:
            self.model_part = self.model.create_model_part(self.output_model_part_name, self.ini_coord[:,0].flatten(),
                                                       self.ini_coord[:,1].flatten(), np.zeros(self.n), np.arange(self.n))

        # interfaces
        self.interface_input = Interface(self.settings['interface_input'], self.model)
        self.interface_input.set_variable_data(self.input_model_part_name, 'heat_flux', self.heat_flux)
        self.interface_output = Interface(self.settings['interface_output'], self.model)
        self.interface_output.set_variable_data(self.output_model_part_name, 'displacement', self.prev_disp + self.dx)

        self.output_solution_step()

    @tools.time_initialize
    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.timestep += 1
        self.prev_disp += self.dx

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # input
        self.interface_input = interface_input.copy()
        self.disp = interface_input.get_variable_data(self.input_model_part_name, 'displacement')
        a = np.pi * (self.d + 2 * self.disp[:, 1]) ** 2 / 4

        # input does not contain boundary conditions
        self.a[1:self.m + 1] = a
        self.a[0] = self.a[1]
        self.a[self.m + 1] = self.a[self.m]

        self.k += 1

        # initial residual
        f = self.get_residual()
        if self.k == 1:
            self.residual0 = np.linalg.norm(f)

        # coupling convergence
        residual = np.linalg.norm(f)
        if self.check_coupling_convergence and residual / self.residual0 < self.newtontol:
            self.coupling_convergence = True
            if self.print_coupling_convergence:
                tools.print_info(f'{self.__class__.__name__} converged')

        # Newton iterations
        converged = False
        if self.residual0:
            for s in range(self.newtonmax):
                j = self.get_jacobian()
                b = -f
                x = solve_banded((self.al, self.au), j, b)
                self.u += x[0::2]
                self.p += x[1::2]
                if self.inlet_variable == 'velocity':
                    self.u[0] = self.get_inlet_boundary()
                elif self.inlet_variable == 'pressure':
                    self.p[0] = self.get_inlet_boundary()
                f = self.get_residual()
                residual = np.linalg.norm(f)
                if residual / self.residual0 < self.newtontol:
                    converged = True
                    break
            if not converged:
                Exception('Newton failed to converge')

        if self.debug:
            file_name = join(self.working_directory, f'input_displacement_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'x-displacement':<22}\t{'y-displacement':<22}"
                           f"\t{'z-displacement':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.disp[i, 0]:<22}\t{self.disp[i, 1]:<22}'
                               f'\t{self.disp[i, 2]:<22}\n')
            p = self.p[1:self.m + 1] * self.rhof
            u = self.u[1:self.m + 1]
            file_name = join(self.working_directory, f'output_pressure_velocity_ts{self.n}_it{self.k}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'velocity':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{p[i]:<22}\t{u[i]:<22}\n')

        # output does not contain boundary conditions
        self.pres = self.p[1:self.m + 1] * self.rhof
        self.interface_output.set_variable_data(self.output_model_part_name, 'pressure', self.pres.reshape(-1, 1))
        self.interface_output.set_variable_data(self.output_model_part_name, 'traction', self.trac)
        return self.interface_output

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

        if self.n > 0 and self.save_restart != 0 and self.n % self.save_restart == 0:
            file_name = join(self.working_directory, f'case_timestep{self.n}.pickle')
            with open(file_name, 'wb') as file:
                pickle.dump({'a': self.a, 'p': self.p, 'u': self.u}, file)
            if self.save_restart < 0 and self.n + self.save_restart > self.timestep_start:
                try:
                    os.remove(join(self.working_directory, f'case_timestep{self.n + self.save_restart}.pickle'))
                except OSError:
                    pass
        if self.debug:
            p = self.p[1:self.m + 1] * self.rhof
            u = self.u[1:self.m + 1]
            file_name = join(self.working_directory, f'pressure_velocity_ts{self.n}.txt')
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'velocity':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{p[i]:<22}\t{u[i]:<22}\n')

    def finalize(self):
        super().finalize()

    def get_inlet_boundary(self):
        if self.inlet_type == 1:
            x = self.inlet_reference \
                + self.inlet_amplitude * np.sin(2 * np.pi * (self.n * self.dt) / self.inlet_period)
        elif self.inlet_type == 2:
            x = self.inlet_reference + (self.inlet_amplitude if self.n <= self.inlet_period / self.dt else 0)
        elif self.inlet_type == 3:
            x = self.inlet_reference \
                + self.inlet_amplitude * (np.sin(np.pi * (self.n * self.dt) / self.inlet_period)) ** 2
        elif self.inlet_type == 4:
            x = self.inlet_reference + self.inlet_amplitude
        else:
            x = self.inlet_reference + self.inlet_amplitude * (self.n * self.dt) / self.inlet_period
        return x

    def get_residual(self):
        usign = self.u[1:self.m + 1] > 0
        ur = self.u[1:self.m + 1] * usign + self.u[2:self.m + 2] * (1.0 - usign)
        ul = self.u[0:self.m] * usign + self.u[1:self.m + 1] * (1.0 - usign)

        f = np.zeros(2 * self.m + 4)
        if self.inlet_variable == 'velocity':
            f[0] = (self.u[0] - self.get_inlet_boundary()) * self.conditioning
            f[1] = (self.p[0] - (2.0 * self.p[1] - self.p[2])) * self.conditioning
        elif self.inlet_variable == 'pressure':
            f[0] = (self.u[0] - (2.0 * self.u[1] - self.u[2])) * self.conditioning
            f[1] = (self.p[0] - self.get_inlet_boundary()) * self.conditioning
        elif self.inlet_variable == 'total_pressure':
            f[0] = (self.u[0] - (2.0 * (self.get_inlet_boundary() - self.p[0])) ** 0.5) * self.conditioning
            f[1] = (self.p[0] - (2.0 * self.p[1] - self.p[2])) * self.conditioning
        f[2:2 * self.m + 2:2] = (self.dz / self.dt * (self.a[1:self.m + 1] - self.an[1:self.m + 1]) * self.unsteady
                                 + (self.u[1:self.m + 1] + self.u[2:self.m + 2])
                                 * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                 - (self.u[1:self.m + 1] + self.u[0:self.m])
                                 * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0
                                 - self.alpha * (self.p[2:self.m + 2] - 2.0 * self.p[1:self.m + 1] + self.p[0:self.m]))
        f[3:2 * self.m + 3:2] = (self.dz / self.dt * (self.u[1:self.m + 1] * self.a[1:self.m + 1]
                                                      - self.un[1:self.m + 1] * self.an[1:self.m + 1]) * self.unsteady
                                 + ur * (self.u[1:self.m + 1] + self.u[2:self.m + 2])
                                 * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                 - ul * (self.u[1:self.m + 1] + self.u[0:self.m])
                                 * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0
                                 + ((self.p[2:self.m + 2] - self.p[1:self.m + 1])
                                    * (self.a[1:self.m + 1] + self.a[2:self.m + 2])
                                    + (self.p[1:self.m + 1] - self.p[0:self.m])
                                    * (self.a[1:self.m + 1] + self.a[0:self.m])) / 4.0)
        f[2 * self.m + 2] = (self.u[self.m + 1] - (2.0 * self.u[self.m] - self.u[self.m - 1])) * self.conditioning
        if self.outlet_type == 1:
            f[2 * self.m + 3] = (self.p[self.m + 1] - (2.0 *
                                                       (self.cmk2 -
                                                        (np.sqrt(self.cmk2 - self.pn[self.m + 1] / 2.0)
                                                         - (self.u[self.m + 1] - self.un[self.m + 1]) / 4.0) ** 2))
                                 ) * self.conditioning
        else:
            f[2 * self.m + 3] = (self.p[self.m + 1] - self.outlet_amplitude) * self.conditioning
        return f

    def get_jacobian(self):
        usign = self.u[1:self.m + 1] > 0
        j = np.zeros((self.al + self.au + 1, 2 * self.m + 4))
        if self.inlet_variable == 'velocity':
            j[self.au + 0 - 0, 0] = 1.0 * self.conditioning  # [0,0]
            j[self.au + 1 - 1, 1] = 1.0 * self.conditioning  # [1,1]
            j[self.au + 1 - 3, 3] = -2.0 * self.conditioning  # [1,3]
            j[self.au + 1 - 5, 5] = 1.0 * self.conditioning  # [1,5]
        elif self.inlet_variable == 'pressure':
            j[self.au + 0 - 0, 0] = 1.0 * self.conditioning  # [0,0]
            j[self.au + 0 - 2, 2] = -2.0 * self.conditioning  # [0,2]
            j[self.au + 0 - 4, 4] = 1.0 * self.conditioning  # [0,4]
            j[self.au + 1 - 1, 1] = 1.0 * self.conditioning  # [1,1]
        elif self.inlet_variable == 'total_pressure':
            j[self.au + 0 - 0, 0] = 1.0 * self.conditioning  # [0,0]
            j[self.au + 0 - 1, 1] = (2.0 * (self.get_inlet_boundary() - self.p[0])) ** -0.5 * self.conditioning  # [1,1]
            j[self.au + 1 - 1, 1] = 1.0 * self.conditioning  # [1,1]
            j[self.au + 1 - 3, 3] = -2.0 * self.conditioning  # [1,3]
            j[self.au + 1 - 5, 5] = 1.0 * self.conditioning  # [1,5]

        j[self.au + 2, 0:2 * self.m + 0:2] = -(self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0  # [2*i, 2*(i-1)]
        j[self.au + 3, 0:2 * self.m + 0:2] = (-((self.u[1:self.m + 1] + 2.0 * self.u[0:self.m]) * usign
                                                + self.u[1:self.m + 1] * (1.0 - usign))
                                              * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i+1, 2*(i-1)]
        j[self.au + 1, 1:2 * self.m + 1:2] = -self.alpha  # [2*i, 2*(i-1)+1]
        j[self.au + 2, 1:2 * self.m + 1:2] = -(self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0  # [2*i+1, 2*(i-1)+1]

        j[self.au + 0, 2:2 * self.m + 2:2] = ((self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                              - (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i, 2*i]
        j[self.au + 1, 2:2 * self.m + 2:2] = (self.dz / self.dt * self.a[1:self.m + 1] * self.unsteady
                                              + ((2.0 * self.u[1:self.m + 1] + self.u[2:self.m + 2]) * usign
                                                 + self.u[2:self.m + 2] * (1.0 - usign))
                                              * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0
                                              - (self.u[0:self.m] * usign
                                                 + (2.0 * self.u[1:self.m + 1] + self.u[0:self.m]) * (1.0 - usign))
                                              * (self.a[1:self.m + 1] + self.a[0:self.m]) / 4.0)  # [2*i+1, 2*i]
        j[self.au - 1, 3:2 * self.m + 3:2] = 2.0 * self.alpha  # [2*i, 2*i+1]
        j[self.au + 0, 3:2 * self.m + 3:2] = (-(self.a[1:self.m + 1] + self.a[2:self.m + 2])
                                              + (self.a[1:self.m + 1] + self.a[0:self.m])) / 4.0  # [2*i+1, 2*i+1]

        j[self.au - 2, 4:2 * self.m + 4:2] = (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0  # [2*i, 2*(i+1)]
        j[self.au - 1, 4:2 * self.m + 4:2] = ((self.u[1:self.m + 1] * usign + (self.u[1:self.m + 1]
                                                                               + 2.0 * self.u[2:self.m + 2])
                                               * (1.0 - usign))
                                              * (self.a[1:self.m + 1] + self.a[2:self.m + 2]) / 4.0)  # [2*i+1, 2*(i+1)]
        j[self.au - 3, 5:2 * self.m + 5:2] = -self.alpha  # [2*i, 2*(i+1)+1]
        j[self.au - 2, 5:2 * self.m + 5:2] = (self.a[1:self.m + 1]
                                              + self.a[2:self.m + 2]) / 4.0  # [2*i+1, 2*(i+1)+1]

        j[self.au + (2 * self.m + 2) - (2 * self.m + 2), 2 * self.m + 2] = 1.0 * self.conditioning  # [2*m+2, 2*m+2]
        j[self.au + (2 * self.m + 2) - (2 * self.m), 2 * self.m] = -2.0 * self.conditioning  # [2*m+2, 2*m]
        j[self.au + (2 * self.m + 2) - (2 * self.m - 2), 2 * self.m - 2] = 1.0 * self.conditioning  # [2*m+2, 2*m-2]
        if self.outlet_type == 1:
            j[self.au + (2 * self.m + 3) - (2 * self.m + 2),
              2 * self.m + 2] = (-(np.sqrt(self.cmk2 - self.pn[self.m + 1] / 2.0)
                                   - (self.u[self.m + 1] - self.un[self.m + 1]) / 4.0)) \
                                * self.conditioning  # [2*m+3, 2*m+2]
        j[self.au + (2 * self.m + 3) - (2 * self.m + 3), 2 * self.m + 3] = 1.0 * self.conditioning  # [2*m+3, 2*m+3]
        return j

    def get_surrogate_jacobian(self):
        # df/d(u,p')
        jf = self.get_jacobian()
        jf = bnd.remove_boundaries(jf, 2)
        jf = jf[1:-1, :]
        jf = bnd.to_dense(jf)

        # dp/df
        jif = np.linalg.inv(jf)
        jif = self.rhof * jif[1::2, :]

        # df/da
        ja = self.get_jacobian_area()

        # da/dr
        jr = 2 * np.sqrt(np.pi * self.a[1: self.m + 1])

        # dp/dr
        jsurrogate = jif @ -ja * jr
        return jsurrogate

    def get_jacobian_area(self):
        usign = self.u[1:self.m + 1] > 0
        ur = self.u[1:self.m + 1] * usign + self.u[2:self.m + 2] * (1.0 - usign)
        ul = self.u[0:self.m] * usign + self.u[1:self.m + 1] * (1.0 - usign)
        j = np.zeros((2 * self.m, self.m))
        for i in range(self.m):
            if i > 0:
                j[2 * i, i - 1] = - (self.u[i] + self.u[i + 1]) / 4  # [2*i, i-1]
                j[2 * i + 1, i - 1] = (- ul[i] * (self.u[i] + self.u[i + 1]) / 4
                                       + (self.p[i + 1] - self.p[i]) / 4)  # [2*i+1, i-1]
            j[2 * i, i] = (self.dz / self.dt * self.unsteady + (self.u[i + 1] + self.u[i + 2]) / 4
                           - (self.u[i] + self.u[i + 1]) / 4)  # [2*i, i]
            j[2 * i + 1, i] = (self.dz / self.dt * self.u[i + 1] * self.unsteady
                               + ur[i] * (self.u[i + 1] + self.u[i + 2]) / 4 - ul[i] * (self.u[i] + self.u[i + 1]) / 4
                               + (self.p[i + 2] - self.p[i + 1]) / 4 + (self.p[i + 1] - self.p[i]) / 4)  # [2*i+1, i]
            if i < self.m - 1:
                j[2 * i, i + 1] = (self.u[i + 1] + self.u[i + 2]) / 4  # [2*i, i+1]
                j[2 * i + 1, i + 1] = (ur[i] * (self.u[i + 1] + self.u[i + 2]) / 4
                                       + (self.p[i + 2] - self.p[i + 1]) / 4)  # [2*i+1, i+1]
        return j
