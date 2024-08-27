from coconut.coupling_components.solver_wrappers.openfoam.openfoam import SolverWrapperOpenFOAM
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

from os.path import join


def create(parameters):
    return SolverWrapperOpenFOAM11(parameters)


class SolverWrapperOpenFOAM11(SolverWrapperOpenFOAM):
    version = '11'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_solver()

        # check that the correct software is available
        self.check_software()

        # raw format
        self.fext = '.xy'  # file extension
        self.nheaderfooter = 1  # number of header and footer lines

        # get solver module
        file_name = join(self.working_directory, 'system/controlDict')
        with open(file_name, 'r') as control_dict_file:
            control_dict = control_dict_file.read()
        self.solver = of_io.get_string(control_dict, 'solver')

    def read_face_centres(self, boundary_name, nfaces):
        filename_x = join(self.working_directory, '0/Ccx')
        filename_y = join(self.working_directory, '0/Ccy')
        filename_z = join(self.working_directory, '0/Ccz')

        x0 = of_io.get_boundary_field(file_name=filename_x, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        y0 = of_io.get_boundary_field(file_name=filename_y, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        z0 = of_io.get_boundary_field(file_name=filename_z, boundary_name=boundary_name, size=nfaces, is_scalar=True)
        return x0, y0, z0

    def wall_shear_stress_dict(self):
        name = self.wall_shear_stress_variable
        dct = (f'    {name}\n'
               f'    {{\n'
               f'        type            {self.wall_shear_stress_variable};\n'
               f'        libs            ("libfieldFunctionObjects.so");\n'
               f'        executeControl  timeStep;\n'
               f'        executeInterval 1;\n'
               f'        writeControl    none;\n'
               f'        patches         $boundaryNames;\n'
               f'    }}\n')
        return dct, name

    def pressure_and_traction_dict(self, boundary_name):
        name = f'coconut_{boundary_name}'
        dct = (f'    {name}\n'
               f'    {{\n'
               f'        type            surfaceFieldValue;\n'
               f'        libs            ("libfieldFunctionObjects.so");\n'
               f'        \n'
               f'        log             false;\n'
               f'        writeControl    timeStep;\n'
               f'        writeInterval   1;\n'
               f'        writeFields     true;\n'
               f'        writeArea       false;\n'
               f'        \n'
               f'        surfaceFormat   raw;\n'
               f'        \n'
               f'        select          patch;\n'
               f'        name            {boundary_name};\n'
               f'        \n'
               f'        operation       none;\n'
               f'        \n'
               f'        fields          (p {self.wall_shear_stress_variable});\n'
               f'    }}\n')
        return dct, name

    def kinematic_conversion(self):
        # based on solver application, set conversion settings from kinematic to static pressure/shear stress
        # typically: incompressible solver, pressure and wallShearStress are kinematic -> multiply with fluid density
        #            compressible solver, pressure and wallShearStress are not kinematic -> do nothing
        kinematic_conversion_dict = {
            'incompressibleFluid': {
                'wall_shear_stress_variable': 'wallShearStress',
                'density_correction_for_pressure': self.settings['density'],
                'density_correction_for_traction': self.settings['density']
            }
        }

        if self.solver not in kinematic_conversion_dict:
            available_applications = ''
            for key in kinematic_conversion_dict:
                available_applications += f'\n\t{key}'
            raise ValueError(f'{self.solver} is not included in the kinematic_conversion_dict '
                             f'used for treatment of (kinematic) pressure and traction\n'
                             f'Add the solver to this dictionary '
                             f'or use one of the existing solvers:{available_applications}')
        else:
            kinematic_conversion = kinematic_conversion_dict[self.solver]
            self.density_for_pressure = kinematic_conversion['density_correction_for_pressure']
            self.density_for_traction = kinematic_conversion['density_correction_for_traction']
            self.wall_shear_stress_variable = kinematic_conversion['wall_shear_stress_variable']
