import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

debug = False

solver_1 = 'fluent.v2023R1'
solver_2 = 'fluent.v2023R1'

# settings
if debug:
    mesh_list = ['M3']
    solid_cores = [2]
    fluid_cores = [4]
    conv_crit_mom = ['1.0E-5']
    conv_crit_energy = ['1.0E-10']
else:
    mesh_list = ['M1', 'M2', 'M3', 'M4']
    solid_cores = [1, 1, 2, 4]
    fluid_cores = [1, 2, 4, 12]
    conv_crit_mom = ['1.0E-5', '1.0E-5', '1.0E-5', '1.0E-5']
    conv_crit_energy = ['1.0E-10', '1.0E-10', '1.0E-10', '1.0E-10']

for i, mesh in enumerate(mesh_list):

    cfd_dir_1 = './' + mesh + '/CFD_1'
    cfd_dir_2 = './' + mesh + '/CFD_2'

    # clean up Fluent
    if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
    if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

    # clean working directories
    shutil.rmtree('./' + mesh, ignore_errors=True)
    # crete new working directory
    os.mkdir('./' + mesh)

    # create new CFD_1 folder
    os.mkdir(cfd_dir_1)
    shutil.copy('setup_files_GC/solid' + "/" + mesh_list[i] + '_solid.msh', cfd_dir_1)
    shutil.copy('setup_files_GC/solid/setup_fluent_1.sh', cfd_dir_1)
    with open('setup_files_GC/solid/case_1.jou') as infile:
        with open(join(cfd_dir_1, 'case_1.jou'), 'w') as outfile:
            for line in infile:
                line = line.replace('|MESH|', mesh)
                line = line.replace('|CC_E|', conv_crit_energy[i])
                outfile.write(line)
    cfd_env_1 = tools.get_solver_env(solver_1, cfd_dir_1)
    subprocess.check_call('./setup_fluent_1.sh', shell=True, cwd=cfd_dir_1, env=cfd_env_1)

    # create new CFD_2 folder
    os.mkdir(cfd_dir_2)
    shutil.copy('setup_files_GC/fluid' + "/" + mesh_list[i] + '_fluid.msh', cfd_dir_2)
    shutil.copy('setup_files_GC/fluid/setup_fluent_2.sh', cfd_dir_2)
    with open('setup_files_GC/fluid/case_2.jou') as infile:
        with open(join(cfd_dir_2, 'case_2.jou'), 'w') as outfile:
            for line in infile:
                line = line.replace('|MESH|', mesh)
                line = line.replace('|CC_M|', conv_crit_mom[i])
                line = line.replace('|CC_E|', conv_crit_energy[i])
                outfile.write(line)
    cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
    subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)

    # copy run_simulation.py script and parameters.json to mesh specific directories
    shutil.copy('setup_files_GC/run_simulation.py', './' + mesh + '/')
    with open('setup_files_GC/parameters.json') as infile:
        with open('./' + mesh + '/parameters.json', 'w') as outfile:
            for line in infile:
                line = line.replace('|MESH|', mesh)
                line = line.replace('|SOLID_CORES|', str(solid_cores[i]))
                line = line.replace('|FLUID_CORES|', str(fluid_cores[i]))
                outfile.write(line)