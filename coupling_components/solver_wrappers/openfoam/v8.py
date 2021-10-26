from coconut.coupling_components.solver_wrappers.openfoam.openfoam import SolverWrapperOpenFOAM
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

from subprocess import check_call
from os.path import join


def create(parameters):
    return SolverWrapperOpenFOAM8(parameters)


class SolverWrapperOpenFOAM8(SolverWrapperOpenFOAM):
    version = '8'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_solver()

        # check that the correct software is available
        self.check_software()

    def write_cell_centres(self):
        check_call('postProcess -func writeCellCentres -time 0 &> log.writeCellCentres;',
                   cwd=self.working_directory, shell=True, env=self.env)

    def read_face_centres(self, boundary_name, nfaces):
        filename_x = join(self.working_directory, '0/Cx')
        filename_y = join(self.working_directory, '0/Cy')
        filename_z = join(self.working_directory, '0/Cz')

        x0 = of_io.get_boundary_field(file_name=filename_x, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        y0 = of_io.get_boundary_field(file_name=filename_y, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        z0 = of_io.get_boundary_field(file_name=filename_z, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        return x0, y0, z0

    def pressure_dict(self, boundary_name):
        dct = (f'PRESSURE_{boundary_name}\n'
               f'{{\n'
               f'type  	             surfaceFieldValue;\n'
               f'libs 	             ("libfieldFunctionObjects.so");\n'
               f'executeControl 	 timeStep;\n'
               f'executeInterval 	 1;\n'
               f'writeControl 	     timeStep;\n'
               f'writeInterval 	     1;\n'
               f'timeFormat 	     fixed;\n'
               f'timePrecision 	     {self.time_precision};\n'
               f'operation 	         none;\n'
               f'writeFields 	     true;\n'
               f'surfaceFormat 	     raw;\n'
               f'regionType 	     patch;\n'
               f'name 	             {boundary_name};\n'
               f'fields              (p);\n'
               f'}}\n')
        return dct

    def wall_shear_stress_dict(self, boundary_name):
        dct = (f'wallShearStress\n'
               f'{{\n'
               f'type  	             wallShearStress;\n'
               f'libs 	             ("libfieldFunctionObjects.so");\n'
               f'executeControl 	 timeStep;\n'
               f'executeInterval 	 1;\n'
               f'writeControl 	     none;\n'
               f'timeFormat          fixed;\n'
               f'timePrecision 	     {self.time_precision};\n'
               f'log 	             false;\n'
               f'}}\n')
        return dct

    def traction_dict(self, boundary_name):
        dct = (f'TRACTION_{boundary_name}\n'
               f'{{\n'
               f'type  	             surfaceFieldValue;\n'
               f'libs 	             ("libfieldFunctionObjects.so");\n'
               f'executeControl 	 timeStep;\n'
               f'executeInterval 	 1;\n'
               f'writeControl 	     timeStep;\n'
               f'writeInterval 	     1;\n'
               f'timeFormat 	     fixed;\n'
               f'timePrecision 	     {self.time_precision};\n'
               f'operation 	         none;\n'
               f'writeFields 	     true;\n'
               f'surfaceFormat 	     raw;\n'
               f'regionType 	     patch;\n'
               f'name 	             {boundary_name};\n'
               f'fields              (wallShearStress);\n'
               f'}}\n')
        return dct
