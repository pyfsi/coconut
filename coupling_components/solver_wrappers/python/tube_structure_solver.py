from coconut.coupling_components.component import Component
from coconut.coupling_components import tools
from coconut.data_structure import Model, Interface

import numpy as np
import os.path as path
from scipy.linalg import solve_banded
import json

def create(parameters):
    return SolverWrapperTubeStructure(parameters)


class SolverWrapperTubeStructure(Component):
    al = 2  # Number of terms below diagonal in matrix
    au = 2  # Number of terms above diagonal in matrix

    def __init__(self, parameters):
        super().__init__()

        # Reading
        self.parameters = parameters
        self.settings = parameters["settings"]
        self.working_directory = self.settings["working_directory"]
        input_file = self.settings["input_file"]
        case_file_name = path.join(self.working_directory, input_file)
        with open(case_file_name, 'r') as case_file:
            new_settings = json.loads(case_file.read())
            self.settings.update(new_settings)  # TODO: inversed priority

        # Settings
        self.unsteady = self.settings.get("unsteady", True)

        l = self.settings["l"]  # Length
        d = self.settings["d"]  # Diameter
        self.rreference = d / 2  # Reference radius of cross section
        self.rhof = self.settings["rhof"]  # Fluid density

        self.preference = self.settings.get("preference", 0)  # Reference pressure

        e = self.settings["e"]  # Young's modulus of structure
        nu = self.settings["nu"]  # Poisson's ratio
        self.h = self.settings["h"]  # Thickness of structure
        self.rhos = self.settings["rhos"]  # Structure density
        self.cmk2 = (e * self.h) / (self.rhof * d)  # Wave speed squared
        self.b1 = (self.h * e) / (1 - nu ** 2) * (self.h ** 2) / 12
        self.b2 = self.b1 * (2 * nu) / self.rreference ** 2
        self.b3 = (self.h * e) / (1 - nu ** 2) * 1 / self.rreference ** 2

        self.m = self.settings["m"]  # Number of segments
        self.dz = l / self.m  # Segment length
        axial_offset = self.settings.get("axial_offset", 0)  # Start position along axis
        self.z = axial_offset + np.arange(self.dz / 2, l, self.dz)  # Data is stored in cell centers

        self.k = 0  # Iteration
        self.n = 0  # Time step (no restart implemented)
        if self.unsteady:
            self.dt = self.settings["delta_t"]  # Time step size
            self.time_discretization = self.settings.get("time_discretization", "backward euler").lower()
            if self.time_discretization == "newmark":
                self.gamma = self.settings["gamma"]  # Newmark parameter: gamma >= 1/2
                self.beta = self.settings["beta"]  # Newmark parameter: beta >= 1/4 * (1/2 + gamma) ^ 2
                self.nm = True
                if not self.gamma >= 0.5 or not self.beta >= 0.25 * (0.5 + self.gamma) ** 2:
                    raise Exception("Inadequate Newmark parameteres")
            elif self.time_discretization == "backward euler":
                self.gamma = 1
                self.beta = 1
                self.nm = False  # used to set rnddot to zero
            else:
                raise ValueError("Time discretization should be 'Newmark' or 'backward Euler'.")
        else:
            # not used
            self.dt = 1  # Time step size default
            self.beta = 1
            self.nm = False

        # Initialization
        self.areference = np.pi * self.rreference ** 2  # Reference area of cross section
        self.p = np.ones(self.m) * self.preference  # Pressure
        self.a = np.ones(self.m) * self.areference  # Area of cross section
        self.r = np.ones(self.m + 4) * self.rreference  # Radius of cross section
        if self.unsteady:
            self.rdot = np.zeros(self.m)  # First derivative of the radius with respect to time in current timestep
            self.rddot = np.zeros(self.m)  # Second derivative of the radius with respect to time in current timestep
        self.rn = np.array(self.r)  # Previous radius of cross section
        self.rndot = np.zeros(self.m)  # First derivative of the radius with respect to time in previous timestep
        self.rnddot = np.zeros(self.m)  # Second derivative of the radius with respect to time in previous timestep

        self.disp = np.zeros((self.m, 3))  # Displacement
        self.trac = np.zeros((self.m, 3))  # Traction (always zero)

        self.conditioning = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.unsteady
                             + 6.0 * self.b1 / self.dz ** 4
                             + 2.0 * self.b2 / self.dz ** 2 + self.b3)  # Factor for conditioning Jacobian

        # ModelParts
        self.model = Model()
        self.model_part = self.model.create_model_part("wall", np.zeros(self.m), self.rreference * np.ones(self.m),
                                                       self.z, np.arange(self.m))  # TODO not use hardcoded mp name

        # Interfaces
        self.interface_input = Interface(self.settings["interface_input"], self.model)
        self.interface_input.set_variable_data("wall", "pressure", self.p)
        self.interface_input.set_variable_data("wall", "traction", self.trac)
        self.interface_output = Interface(self.model, self.settings["interface_output"])
        self.interface_input.set_variable_data("wall", "displacement", self.disp)

        # run time
        self.run_time = 0.0

        # Debug
        self.debug = False  # Set on True to save input and output of each iteration of every time step
        self.output_solution_step()

    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        self.k = 0
        self.n += 1
        if self.unsteady:
            self.rn = np.array(self.r)
            self.rndot = np.array(self.rdot)
            self.rnddot = np.array(self.rddot)

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        # Input
        self.interface_input = interface_input.copy()
        self.p = interface_input.get_variable_data("wall", "pressure")

        # Solve system
        f = self.get_residual()
        residual0 = np.linalg.norm(f)
        if residual0:
            j = self.get_jacobian()
            b = -f
            x = solve_banded((self.al, self.au), j, b)
            self.r += x

        self.k += 1
        if self.debug:
            file_name = self.working_directory + f"/Input_Pressure_Traction_TS{self.n}_IT{self.k}.txt"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'pressure':<22}\t{'x-traction':<22}"
                           f"\t{'y-traction':<22}\t{'z-traction':<22}\n")
                for i in range(len(self.z)):
                    file.write(f"{self.z[i]:<22}\t{self.p[i]:<22}\t{self.trac[i, 0]:<22}"
                               f"\t{self.trac[i, 1]:<22}\t{self.trac[i, 2]:<22}\n")
            file_name = self.working_directory + f"/output_Area_TS{self.n}_IT{self.k}.txt"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

        # output does not contain boundary conditions
        self.a = self.r[2:self.m + 2] ** 2 * np.pi
        self.disp[:, 1] = self.r[2:self.m + 2] - self.rreference
        self.interface_output.set_variable_data("wall", "pressure", self.disp)
        return self.interface_output  # TODO: make copy?

    def finalize_solution_step(self):
        super().finalize_solution_step()

        if self.unsteady:
            self.rddot = ((self.r[2:self.m + 2] - self.rn[2:self.m + 2]) / (self.beta * self.dt ** 2)
                          - self.rndot / (self.beta * self.dt) - self.rnddot * (1 / (2 * self.beta) - 1) * self.nm)
            self.rdot = self.rndot + self.dt * (1 - self.gamma) * self.rnddot + self.dt * self.gamma * self.rddot

    def finalize(self):
        super().finalize()

    def output_solution_step(self):
        if self.debug:
            file_name = self.working_directory + f"/Area_TS{self.n}.txt"
            with open(file_name, 'w') as file:
                file.write(f"{'z-coordinate':<22}\t{'area':<22}\n")
                for i in range(len(self.z)):
                    file.write(f'{self.z[i]:<22}\t{self.a[i]:<22}\n')

    def get_interface_input(self):  # TODO: need to have latest data?
        return self.interface_input.copy()

    def set_interface_input(self):  # TODO: remove?
        raise Exception("This solver interface provides no mapping.")

    def get_interface_output(self):  # TODO: need to have latest data?
        return self.interface_output.copy()

    def set_interface_output(self):  # TODO: remove?
        raise Exception("This solver interface provides no mapping.")

    def get_residual(self):
        f = np.zeros(self.m + 4)
        f[0] = (self.r[0] - self.rreference) * self.conditioning
        f[1] = (self.r[1] - self.rreference) * self.conditioning
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
                           * self.unsteady)
        f[self.m + 2] = (self.r[self.m + 2] - self.rreference) * self.conditioning
        f[self.m + 3] = (self.r[self.m + 3] - self.rreference) * self.conditioning
        return f

    def get_jacobian(self):
        j = np.zeros((self.al + self.au + 1, self.m + 4))
        j[self.au + 0 - 0, 0] = 1.0 * self.conditioning  # [0, 0]
        j[self.au + 1 - 1, 1] = 1.0 * self.conditioning  # [1, 1]
        j[self.au + 2, 0:self.m] = self.b1 / self.dz ** 4  # [i, (i - 2)]
        j[self.au + 1, 1:self.m + 1] = - 4.0 * self.b1 / self.dz ** 4 - self.b2 / self.dz ** 2  # [i, (i - 1)]
        j[self.au + 0, 2:self.m + 2] = ((self.rhos * self.h) / (self.beta * self.dt ** 2) * self.unsteady
                                        + 6.0 * self.b1 / self.dz ** 4 + 2.0 * self.b2 / self.dz ** 2
                                        + self.b3)  # [i, i]
        j[self.au - 1, 3:self.m + 3] = - 4.0 * self.b1 / self.dz ** 4 - self.b2 / self.dz ** 2  # [i, (i + 1)]
        j[self.au - 2, 4:self.m + 4] = self.b1 / self.dz ** 4  # [i, (i + 2)]
        j[self.au + (self.m + 2) - (self.m + 2), self.m + 2] = 1.0 * self.conditioning  # [m + 2, m + 2]
        j[self.au + (self.m + 3) - (self.m + 3), self.m + 3] = 1.0 * self.conditioning  # [m + 3, m + 3]
        return j
