from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam_base.openfoam_base import SolverWrapperOpenFoamBase

import os
import subprocess


def create(parameters):
    return SolverWrapperOpenFoamORG(parameters)


class SolverWrapperOpenFoamORG(SolverWrapperOpenFoamBase):

    @tools.time_initialize
    def __init__(self, parameters):
        super().__init__(parameters)

    @tools.time_initialize
    def initialize(self):
        super().initialize()

    def initialize_solution_step(self):
        super().initialize_solution_step()

    @tools.time_solve_solution_step
    def solve_solution_step(self, interface_input):
        out = super().solve_solution_step(interface_input)
        return out

    def finalize_solution_step(self):
        super().finalize_solution_step()

    @tools.time_save
    def output_solution_step(self):
        super().output_solution_step()

    def finalize(self):
        super().finalize()

    def compile_adapted_openfoam_solver(self):
        # compile openfoam adapted solver
        solver_dir = os.path.join(os.path.dirname(__file__), f'v{self.version.replace(".", "")}', self.application)
        try:
            if self.compile_clean:
                subprocess.check_call(f'wclean {solver_dir} && wmake {solver_dir} &> log.wmake',
                                      cwd=self.working_directory, shell=True,
                                      env=self.env)
            else:
                subprocess.check_call(f'wmake {solver_dir} &> log.wmake', cwd=self.working_directory, shell=True,
                                      env=self.env)
        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'Compilation of {self.application} failed. Check {os.path.join(self.working_directory, "log.wmake")}')

    def check_software(self):
        try:
            version_nr = subprocess.check_output('echo $WM_PROJECT_VERSION',
                                                 shell=True, env=self.env).decode('utf-8').strip()  # output b'XX\n'
        except subprocess.CalledProcessError:
            raise RuntimeError(
                f'OpenFOAM not loaded properly. Check if the solver load commands for the "machine_name"'
                f' are correct in solver_modules.py.')

        # check version
        if version_nr != self.version:
            raise RuntimeError(
                f'OpenFOAM-{self.version} should be loaded! Currently, OpenFOAM-{version_nr} is loaded.'
                f' Check if the solver load commands for the "machine_name" are correct in solver_modules.py.')
