import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

debug = False

solver_1 = 'cht_fluent.cht_v2023R1'
solver_2 = 'cht_fluent.cht_v2023R1'

# settings
mesh = 'M3'
solid_cores = 1
fluid_cores = 4
conv_crit_mom = '1.0E-5'
conv_crit_energy = '1.0E-10'
conduction = 0.2876 # W/mK

# relaxation factors
omega_gs = [[0.4, 0.65, 0.75, 0.75], [0.48, 0.32, 0.22, 0.20]]
omega_aitken = [[0.5, 0.65, 0.85, 0.95], [0.6, 0.5, 0.4, 0.4]]

def relax(v, i, method):
    if method == 'iqni':
        return 1.0
    elif method == 'relaxation':
        return omega_gs[i][v]
    elif method == 'aitken':
        return omega_aitken[i][v]
    else:
        raise ValueError(f'Method {method} not implemented.')

if debug:
    velocity_list = [24] # m/s
    velocity_names = ['24m_s']
    solver_order = ['FFTB']
    methods = ['iqni']
else:
    velocity_list = [2.4, 12, 24, 36] # m/s
    velocity_names = ['2-4m_s', '12m_s', '24m_s', '36m_s']
    solver_order = ['TFFB', 'FFTB']
    methods = ['relaxation', 'aitken', 'iqni']

for v, velocity in enumerate(velocity_list):
    for i, order in enumerate(solver_order):
        for k, method in enumerate(methods):

            cfd_dir_1 = './' + velocity_names[v] + '/' + order + '_' + method + '/CFD_1'
            cfd_dir_2 = './' + velocity_names[v] + '/' + order + '_' + method + '/CFD_2'

            # clean up Fluent
            if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
                subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
            if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
                subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

for v, velocity in enumerate(velocity_list):
    # clean working directories and create new ones
    shutil.rmtree('./' + velocity_names[v], ignore_errors=True)
    os.mkdir('./' + velocity_names[v])

    for i, order in enumerate(solver_order):
        for k, method in enumerate(methods):
            os.mkdir('./' + velocity_names[v] + '/' + order + '_' + method)

            cfd_dir_1 = './' + velocity_names[v] + '/' + order + '_' + method + '/CFD_1'
            cfd_dir_2 = './' + velocity_names[v] + '/' + order + '_' + method + '/CFD_2'

            # create new CFD_1 folder
            os.mkdir(cfd_dir_1)
            shutil.copy('setup_files/solid' + "/" + mesh + '_solid.msh', cfd_dir_1)
            shutil.copy('setup_files/solid/setup_fluent_1.sh', cfd_dir_1)
            with open('setup_files/solid/case_1.jou') as infile:
                with open(join(cfd_dir_1, 'case_1.jou'), 'w') as outfile:
                    for line in infile:
                        line = line.replace('|MESH|', mesh)
                        line = line.replace('|CC_E|', conv_crit_energy)
                        line = line.replace('|COND|', str(conduction))
                        outfile.write(line)
            cfd_env_1 = tools.get_solver_env(solver_1, cfd_dir_1)
            subprocess.check_call('./setup_fluent_1.sh', shell=True, cwd=cfd_dir_1, env=cfd_env_1)

            # create new CFD_2 folder
            os.mkdir(cfd_dir_2)
            shutil.copy('setup_files/fluid' + "/" + mesh + '_fluid.msh', cfd_dir_2)
            shutil.copy('setup_files/fluid/setup_fluent_2.sh', cfd_dir_2)
            with open('setup_files/fluid/case_2.jou') as infile:
                with open(join(cfd_dir_2, 'case_2.jou'), 'w') as outfile:
                    for line in infile:
                        line = line.replace('|MESH|', mesh)
                        line = line.replace('|CC_M|', conv_crit_mom)
                        line = line.replace('|CC_E|', conv_crit_energy)
                        line = line.replace('|VEL|', str(velocity))
                        outfile.write(line)
            cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
            subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)

            # copy run_simulation.py script and parameters.json to mesh specific directories
            shutil.copy('setup_files/run_simulation.py', './' + velocity_names[v] + '/' + order + '_' + method + '/')
            with open('setup_files/parameters_' + order + '.json') as infile:
                with open('./' + velocity_names[v] + '/' + order + '_' + method + '/parameters.json', 'w') as outfile:
                    for line in infile:
                        line = line.replace('|MESH|', mesh)
                        line = line.replace('|SOLID_CORES|', str(solid_cores))
                        line = line.replace('|FLUID_CORES|', str(fluid_cores))
                        line = line.replace('|METHOD|', method)
                        line = line.replace('|OMEGA|', str(relax(v, i, method)))
                        line = line.replace('|OMEGA_MAX|', str(relax(v, i, method)))
                        outfile.write(line)