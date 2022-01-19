from coconut.coupling_components.solver_wrappers.openfoam.openfoam import SolverWrapperOpenFOAM
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

from subprocess import check_call
from os.path import join


def create(parameters):
    return SolverWrapperOpenFOAM41(parameters)


class SolverWrapperOpenFOAM41(SolverWrapperOpenFOAM):
    version = '4.1'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_solver()

        # check that the correct software is available
        self.check_software()

    def write_cell_centres(self):
        check_call('writeCellCentres -time 0 &> log.writeCellCentres;', cwd=self.working_directory, shell=True,
                   env=self.env)

    def read_face_centres(self, boundary_name, nfaces):
        filename_x = join(self.working_directory, '0/ccx')
        filename_y = join(self.working_directory, '0/ccy')
        filename_z = join(self.working_directory, '0/ccz')

        x0 = of_io.get_boundary_field(file_name=filename_x, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        y0 = of_io.get_boundary_field(file_name=filename_y, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        z0 = of_io.get_boundary_field(file_name=filename_z, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        return x0, y0, z0

    def pressure_dict(self, boundary_name, pre=''):
        dct = (f'{pre}PRESSURE_{boundary_name}\n'
               f'{pre}{{\n'
               f'{pre}    type            surfaceRegion;\n'
               f'{pre}    libs            ("libfieldFunctionObjects.so");\n'
               f'{pre}    executeControl  timeStep;\n'
               f'{pre}    executeInterval 1;\n'
               f'{pre}    writeControl    timeStep;\n'
               f'{pre}    writeInterval   1;\n'
               f'{pre}    timeFormat      fixed;\n'
               f'{pre}    timePrecision   {self.time_precision};\n'
               f'{pre}    operation       none;\n'
               f'{pre}    writeFields     true;\n'
               f'{pre}    surfaceFormat   raw;\n'
               f'{pre}    regionType      patch;\n'
               f'{pre}    name            {boundary_name};\n'
               f'{pre}    fields          (p);\n'
               f'{pre}}}\n')
        return dct

    def wall_shear_stress_dict(self, boundary_name, pre=''):
        dct = (f'{pre}wallShearStress\n'
               f'{pre}{{\n'
               f'{pre}    type            wallShearStress;\n'
               f'{pre}    libs            ("libfieldFunctionObjects.so");\n'
               f'{pre}    executeControl  timeStep;\n'
               f'{pre}    executeInterval 1;\n'
               f'{pre}    writeControl    none;\n'
               f'{pre}    timeFormat      fixed;\n'
               f'{pre}    timePrecision   {self.time_precision};\n'
               f'{pre}    log             false;\n'
               f'{pre}}}\n')
        return dct

    def traction_dict(self, boundary_name, pre=''):
        dct = (f'{pre}TRACTION_{boundary_name}\n'
               f'{pre}{{\n'
               f'{pre}    type            surfaceRegion;\n'
               f'{pre}    libs            ("libfieldFunctionObjects.so");\n'
               f'{pre}    executeControl  timeStep;\n'
               f'{pre}    executeInterval 1;\n'
               f'{pre}    writeControl    timeStep;\n'
               f'{pre}    writeInterval   1;\n'
               f'{pre}    timeFormat      fixed;\n'
               f'{pre}    timePrecision   {self.time_precision};\n'
               f'{pre}    operation       none;\n'
               f'{pre}    writeFields     true;\n'
               f'{pre}    surfaceFormat   raw;\n'
               f'{pre}    regionType      patch;\n'
               f'{pre}    name            {boundary_name};\n'
               f'{pre}    fields          (wallShearStress);\n'
               f'{pre}}}\n')
        return dct
