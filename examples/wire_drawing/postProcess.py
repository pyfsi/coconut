import numpy as np
from coconut.coupling_components.solver_wrappers.openfoam import openfoam_io as of_io
from coconut import data_structure
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import host_subplot
# import mpl_toolkits.axisartist as AA
from scipy import interpolate
import pandas as pd
import json
import os
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

def create(parameters):
    return postProcess(parameters)

class postProcess():

    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters['coupled_solver']
        self.wrapper = self.settings['solver_wrappers']
        self.structural_wrapper = self.wrapper[1]
        self.settings_structure = self.structural_wrapper['settings']
        self.wrapper = self.settings_structure['solver_wrapper']
        self.settings_wrapper = self.wrapper['settings']

        self.die_min = self.settings_wrapper['axial_coordinate_die_min']
        self.die_max = self.settings_wrapper['axial_coordinate_die_max']
        self.time_folder = self.settings_wrapper['time_folder_postProcess']

        cwd = os.getcwd()

        CSM_directory = cwd + '/CSM'
        CFD_directory = cwd + '/CFD'

        boundary_name_CFD = 'dieWall'
        boundary_name_CSM = 'wireTopContact'

        node_ids_CFD, node_coords_CFD = of_io.get_boundary_points(case_directory=CFD_directory, time_folder=self.time_folder,
                                                                      boundary_name=boundary_name_CFD)

        color1 = 'r'
        color2 = 'b'
        color3 = 'g'
        color4 = 'c'
        color5 = plt.cm.viridis(0)

        # postprocess pressure and traction distribution
        pressure=[]
        shearstress=[]

        #pressure
        file_name_pressure = os.path.join('CFD/postProcessing/PRESSURE_Wire/surface/' + self.time_folder,
                                          f'p_patch_Wire.raw')
        pressureList = np.genfromtxt(file_name_pressure)
        pressure.append(pressureList)

        pressureList = np.genfromtxt (file_name_pressure)
        pressure.append(pressureList)

        #shear stress
        file_name_shear = os.path.join('CFD/postProcessing/TRACTION_Wire/surface/' + self.time_folder,
                                          f'wallShearStress_patch_Wire.raw')
        shearStressList = np.genfromtxt(file_name_shear)
        shearstress.append(shearStressList)
        wall_shear_stress = shearstress[0]
        wall_shear_stress[:,3] = np.sqrt(wall_shear_stress[:,3]**2 + wall_shear_stress[:,4]**2)

        # plot
        plt.style.use('default')

        fig, axb1 = plt.subplots(figsize=(7, 4.5))
        axb2 = axb1.twinx()
        axb3 = axb1.twinx()
        axb2.ticklabel_format(style='plain', axis='y', scilimits=None)
        axb3.ticklabel_format(style='plain', axis='y')

        axb1.set_ylim(0.4, 0.9)
        axb2.set_ylim(0, 1.8)
        axb3.set_ylim(0, 0.1)

        start, end = axb1.get_ylim()
        axb1.yaxis.set_ticks(np.arange(start, end, 0.05))

        axb1.set_xlabel('Axial wire-die position (mm)', fontsize=12)
        axb1.set_ylabel('Radial position (mm)', fontsize=12)
        axb2.set_ylabel('Pressure (GPa)', fontsize=12)
        axb3.set_ylabel('Shear stress (GPa)', fontsize=12)

        axb1.tick_params(axis='x', labelsize=12)
        axb1.tick_params(axis='y', labelsize=12)
        axb2.tick_params(axis='y', labelsize=12)
        axb3.tick_params(axis='y', labelsize=12)
        axb1.xaxis.set_minor_locator(AutoMinorLocator(5))
        axb1.yaxis.set_minor_locator(AutoMinorLocator(4))
        axb2.yaxis.set_minor_locator(AutoMinorLocator(5))
        axb3.yaxis.set_minor_locator(AutoMinorLocator(5))

        p1, = axb1.plot(node_coords_CFD[:, 0] * 10e2, node_coords_CFD[:, 1] * 10e2, color=color1, label="Die position")
        pressure_plot = pressure[0]
        p2, = axb2.plot(pressure_plot[:, 0] * 10e2, pressure_plot[:, 3] / (1e9), color=color4,
                        label=f'Pressure distribution soap')
        p3, = axb3.plot(wall_shear_stress[:, 0] * 10e2, wall_shear_stress[:, 3] / (1e9), color=color5,
                        label='Shear stress distribution soap')

        lns = [p1, p2, p3]
        axb1.legend(handles=lns, loc= 6)

        axb3.spines['right'].set_position(('outward', 50))

        axb1.yaxis.label.set_color(p1.get_color())
        axb2.yaxis.label.set_color(p2.get_color())
        axb3.yaxis.label.set_color(p3.get_color())

        fig.tight_layout()

        # Yield and Cauchy Stress

        # specify location of Yield and Cauchy Stress
        self.model = data_structure.Model()
        mp_name = f'{boundary_name_CSM}_output'
        node_ids, node_coords = of_io.get_boundary_points(case_directory=CSM_directory,
                                                          time_folder=self.time_folder,
                                                          boundary_name=boundary_name_CSM)

        self.model.create_model_part(mp_name, node_coords[:, 0], node_coords[:, 1],
                                     node_coords[:, 2], node_ids)

        mp = self.model.get_model_part(mp_name)
        # TODO: Be carefull for nfaces, it is not the correct application compared with solver-wrappers!!!
        nfaces = int(mp.size / 2 - 1)

        filename_sigmaCauchy = os.path.join(CSM_directory, self.time_folder, 'sigmaCauchyEq')
        filename_sigmaY = os.path.join(CSM_directory, self.time_folder, 'sigmaY')
        filename_rad_disp = os.path.join(CSM_directory, self.time_folder, 'U')


        cauchy_field = of_io.get_boundary_field(file_name=filename_sigmaCauchy, boundary_name=boundary_name_CSM,
                                                size=nfaces, is_scalar=True)
        sigmaY_field = of_io.get_boundary_field(file_name=filename_sigmaY, boundary_name=boundary_name_CSM,
                                                size=nfaces, is_scalar=True)
        rad_disp = of_io.get_boundary_field(file_name=filename_rad_disp, boundary_name=boundary_name_CSM,
                                                size=nfaces, is_scalar=False)

        filename_ccx = os.path.join(CSM_directory, self.time_folder, 'ccx')
        filename_ccy = os.path.join(CSM_directory, self.time_folder, 'ccy')
        filename_ccz = os.path.join(CSM_directory, self.time_folder, 'ccz')
        x = of_io.get_boundary_field(file_name=filename_ccx, boundary_name=boundary_name_CSM, size=nfaces,
                                     is_scalar=True)
        y = of_io.get_boundary_field(file_name=filename_ccy, boundary_name=boundary_name_CSM, size=nfaces,
                                     is_scalar=True)
        z = of_io.get_boundary_field(file_name=filename_ccz, boundary_name=boundary_name_CSM, size=nfaces,
                                     is_scalar=True)

        arr_test = np.array([x, cauchy_field, sigmaY_field, rad_disp[:,1],y])
        sorted = arr_test[:, arr_test[0].argsort()]


        # plot
        plt.style.use('default')

        fig, axc1 = plt.subplots(figsize=(6.125, 4.5))
        axc2 = axc1.twinx()
        axc1.ticklabel_format(style='sci', axis='y', scilimits=None)
        axc2.ticklabel_format(style='sci', axis='y', scilimits=None)

        axc1.set_xlim(-3, 3)
        axc1.set_ylim(0.4, 0.9)
        axc2.set_ylim(0, 2)

        axc1.set_xlabel('Axial wire-die position (mm)', fontsize=12)
        axc1.set_ylabel('Radial position (mm)', fontsize=12)
        axc2.set_ylabel('Stress (GPa)', fontsize=12)

        axc1.tick_params(axis='x', labelsize=12)
        axc1.tick_params(axis='y', labelsize=12)
        axc2.tick_params(axis='y', labelsize=12)

        start, end = axc1.get_ylim()
        axc1.yaxis.set_ticks(np.arange(start, end, 0.05))

        axc1.xaxis.set_minor_locator(AutoMinorLocator(5))
        axc1.yaxis.set_minor_locator(AutoMinorLocator(4))
        axc2.yaxis.set_minor_locator(AutoMinorLocator(5))

        color6 = plt.cm.viridis(0.8)
        color7 = "k"

        p1, = axc1.plot(node_coords_CFD[:, 0] * 10e2, node_coords_CFD[:, 1] * 10e2, color=color1, label="Die position")
        p2, = axc2.plot(sorted[0] * 10e2, sorted[1] / (10e8), color=color6, label=f'Equivalent Cauchy stress at top wire')
        p3, = axc2.plot(sorted[0] * 10e2, sorted[2]/10e8 , color=color7, linestyle='dashed',
                        label='Yield stress at top wire')

        lns = [p1, p2, p3]
        axc1.legend(handles=lns, loc=1)

        fig.tight_layout()

        fig.savefig(f'{self.time_folder}_stress.png', dpi=800)

        # plot radial_disp and film thickness oter apporach to avoid problems with layering
        node_ids_CFD, node_coords_CFD = of_io.get_boundary_points(case_directory=CFD_directory,
                                                                  time_folder=self.time_folder,
                                                                  boundary_name=boundary_name_CFD)
        q = int(node_coords_CFD[:, 0].size / 2)
        filter_filter_node_coords_CFD = np.zeros([q, 3])
        m = 0
        for i in range(node_coords_CFD[:, 0].size):
            if node_coords_CFD[i, 2] > 0:
                filter_filter_node_coords_CFD[m, 0] = node_coords_CFD[i, 0]
                filter_filter_node_coords_CFD[m, 1] = node_coords_CFD[i, 1]
                filter_filter_node_coords_CFD[m, 2] = node_coords_CFD[i, 2]
                m += 1

        f = interpolate.interp1d(sorted[0], sorted[4])

        xnew = filter_filter_node_coords_CFD[:, 0]
        ynew = f(xnew)
        diff = filter_filter_node_coords_CFD[:,1] - ynew

        plt.style.use('default')

        fig, ax1 = plt.subplots(figsize=(8,5))
        ax2 = ax1.twinx()
        ax3 = ax1.twinx()

        ax1.set_xlim(-1, 1.5)
        ax1.set_ylim(0.6, 0.8)
        ax2.set_ylim(-10, 60)
        ax3.set_ylim(-50, 125)

        ax1.set_xlabel('Axial wire-die position (mm)', fontsize=12)
        ax1.set_ylabel('Radial position (mm)', fontsize=12)
        ax2.set_ylabel('film thickness (µm)', fontsize=12)
        ax3.set_ylabel('Displacement (µm)', fontsize=12)

        ax1.tick_params(axis='x', labelsize=12)
        ax1.tick_params(axis='y', labelsize=12)
        ax2.tick_params(axis='y', labelsize=12)
        ax3.tick_params(axis='y', labelsize=12)
        ax1.xaxis.set_minor_locator(AutoMinorLocator(4))
        ax1.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
        ax3.yaxis.set_minor_locator(AutoMinorLocator(5))

        p1, = ax1.plot(node_coords_CFD[:, 0]*10e2 , node_coords_CFD[:, 1]*10e2, color=color1, label="Die position")
        p2, = ax2.plot(xnew*10e2, diff*10e5, color=color2,
                        label=f'Film thickness soap')
        p3, = ax3.plot(sorted[0]*10e2, sorted[3]*10e5, color=color3,
                        label='Radial displacement wire')

        lns = [p1, p2,p3]
        ax1.legend(handles=lns, loc= 'best')
        ax3.spines['right'].set_position(('outward', 50))

        ax1.yaxis.label.set_color(p1.get_color())
        ax2.yaxis.label.set_color(p2.get_color())
        ax3.yaxis.label.set_color(p3.get_color())

        fig.tight_layout()

        plt.show()

# Import parameters
parameter_file_name = "parameters.json"
with open(parameter_file_name, 'r') as parameter_file:
   parameters = json.load(parameter_file)

A = postProcess(parameters)

