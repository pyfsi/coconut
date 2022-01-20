from coconut.coupling_components.solver_wrappers.openfoam.openfoam import SolverWrapperOpenFOAM
from coconut import tools
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io

from subprocess import check_call, CalledProcessError
from os.path import join
import numpy as np
import os


def create(parameters):
    return SolverWrapperOpenFOAM41Water(parameters)


class SolverWrapperOpenFOAM41Water(SolverWrapperOpenFOAM):
    version = '4.1'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.application = 'coconut_interFoam'
        self.env = tools.get_solver_env("coconut.coupling_components.solver_wrappers.openfoam.v41", self.working_directory)

        # compile adapted openfoam software
        self.compile_adapted_openfoam_solver()

        # check that the correct software is available
        self.check_software()

    def compile_adapted_openfoam_solver(self):
        # compile openfoam adapted solver
        solver_dir = os.path.join(os.path.dirname(__file__), 'v41', self.application)
        try:
            check_call(f'wmake {solver_dir} &> log.wmake', cwd=self.working_directory, shell=True,
                                  env=self.env)
        except CalledProcessError:
            raise RuntimeError(
                f'Compilation of {self.application} failed. Check {os.path.join(self.working_directory, "log.wmake")}')

    def write_cell_centres(self):
        check_call('writeCellCentres -time 0 &> log.writeCellCentres;', cwd=self.working_directory, shell=True,
                   env=self.env)

    def read_face_centres(self, boundary_name, nfaces):
        filename_x = join(self.working_directory, '0/ccx')
        filename_y = join(self.working_directory, '0/ccy')
        filename_z = join(self.working_directory, '0/ccz')

        x0 = of_io.get_boundary_field(file_name=filename_x, boundary_name=boundary_name, size=nfaces,
                                      is_scalar=True)
        y0 = of_io.get_boundary_field(file_name=filename_y, boundary_name=boundary_name, size=nfaces,
                                      is_scalar=True)
        z0 = of_io.get_boundary_field(file_name=filename_z, boundary_name=boundary_name, size=nfaces,
                                      is_scalar=True)
        return x0, y0, z0

    def pressure_dict(self, boundary_name):
        dct = (f'PRESSURE_{boundary_name}\n'
               f'{{\n'
               f'type  	             surfaceRegion;\n'
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
               f'fields              (p alpha.water);\n'
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
               f'type  	             surfaceRegion;\n'
               f'libs 	             ("libfieldFunctionObjects.so");\n'
               f'executeControl      timeStep;\n'
               f'executeInterval 	 1;\n'
               f'writeControl 	     timeStep;\n'
               f'writeInterval 	     1;\n'
               f'timeFormat 	     fixed;\n'
               f'timePrecision 	     {self.time_precision};\n'
               f'operation 	         none;\n'
               f'writeFields 	     true;\n'
               f'surfaceFormat       raw;\n'
               f'regionType 	     patch;\n'
               f'name 	             {boundary_name};\n'
               f'fields              (wallShearStress);\n'
               f'}}\n')
        return dct

    def read_node_output(self):
        """
        reads the pressure, alpha and wall shear stress from the <case directory>/postProcessing for serial and parallel. In
        case of parallel, it uses mp.sequence (using faceProcAddressing) to map the values to the face centres.

        :return:
        """
        alpha_thershold = 0.5
        for boundary in self.boundary_names:
            # specify location of pressure and traction
            traction_name = 'TRACTION_' + boundary
            pressure_name = 'PRESSURE_' + boundary
            mp_name = f'{boundary}_output'
            mp = self.model.get_model_part(mp_name)
            nfaces = mp.size
            wss_filename = os.path.join(self.working_directory, 'postProcessing', traction_name, 'surface',
                                        self.cur_timestamp, 'wallShearStress_patch_' + boundary + '.raw')
            pres_filepath = os.path.join(self.working_directory, 'postProcessing', pressure_name, 'surface',
                                         self.cur_timestamp, 'p_patch_' + boundary + '.raw')
            alpha_filepath = os.path.join(self.working_directory, 'postProcessing', pressure_name, 'surface',
                                         self.cur_timestamp, 'alpha.water_patch_' + boundary + '.raw')

            # check if the pressure and wall shear stress files completed by openfoam and read data
            self.check_output_file(wss_filename, nfaces)
            wss_tmp = np.loadtxt(wss_filename, comments='#')[:, 3:]
            self.check_output_file(pres_filepath, nfaces)
            pres_tmp = np.loadtxt(pres_filepath, comments='#')[:, 3]
            self.check_output_file(alpha_filepath, nfaces)
            alpha_tmp = np.loadtxt(alpha_filepath, comments='#')[:, 3]
            #alpha_tmp = np.where(alpha_tmp<0.1, 0.0, 1.0)

            if self.settings['parallel']:
                pos_list = self.mp_out_reconstruct_seq_dict[mp_name]
            else:
                pos_list = [pos for pos in range(0, nfaces)]

            wall_shear_stress = np.empty_like(wss_tmp)
            pressure = np.empty((pres_tmp.size, 1))
            alpha = np.empty((alpha_tmp.size, 1))

            wall_shear_stress[pos_list,] = wss_tmp[:, ]
            pressure[pos_list, 0] = pres_tmp
            alpha[pos_list, 0] = alpha_tmp

            self.interface_output.set_variable_data(mp_name, 'traction', wall_shear_stress * -1 * alpha)
            self.interface_output.set_variable_data(mp_name, 'pressure', pressure * alpha)

        def delete_prev_iter_output(self):
            # pressure and wall shear stress files are removed to avoid openfoam to append data in the new iteration
            for boundary in self.boundary_names:
                # specify location of pressure and traction
                traction_name = 'TRACTION_' + boundary
                pressure_name = 'PRESSURE_' + boundary
                wss_file = os.path.join(self.working_directory, 'postProcessing', traction_name, 'surface',
                                        self.cur_timestamp, 'wallShearStress_patch_' + boundary + '.raw')
                pres_file = os.path.join(self.working_directory, 'postProcessing', pressure_name, 'surface',
                                         self.cur_timestamp, 'p_patch_' + boundary + '.raw')
                alpha_file = os.path.join(self.working_directory, 'postProcessing', pressure_name, 'surface',
                                         self.cur_timestamp, 'alpha.water_patch_' + boundary + '.raw')
                if os.path.isfile(wss_file):
                    os.remove(wss_file)
                if os.path.isfile(pres_file):
                    os.remove(pres_file)
                if os.path.isfile(alpha_file):
                    os.remove(alpha_file)
