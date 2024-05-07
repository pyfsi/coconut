import glob
import shutil
import subprocess
from os.path import join
from os import remove

from coconut import tools

solver_1 = 'fluent.v2023R1'
solver_2 = 'fluent.v2023R1'
mesh_list = ['M1', 'M2', 'M3', 'M4']

for i, mesh in enumerate(mesh_list):

    cfd_dir_1 = './' + mesh + '/CFD_1'
    cfd_dir_2 = './' + mesh + '/CFD_2'

    # clean up Fluent
    if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
    if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

    # clean working directories
    shutil.rmtree(cfd_dir_1, ignore_errors=True)
    shutil.rmtree(cfd_dir_2, ignore_errors=True)

    # create new CFD_1 folder
    shutil.copytree('setup_files_GC/solid', cfd_dir_1)
    for j in range(len(mesh_list)):
        if j != i: remove(cfd_dir_1 + "/" + mesh_list[j] + '_solid.msh')
    with open(join(cfd_dir_1, 'case_1.jou')) as infile:
        with open(join(cfd_dir_1, 'Case_1.jou'), 'w') as outfile:
            for line in infile:
                line = line.replace('|MESH|', mesh)
                outfile.write(line)
    remove(join(cfd_dir_1, 'case_1.jou'))
    cfd_env_1 = tools.get_solver_env(solver_1, cfd_dir_1)
    subprocess.check_call('./setup_fluent_1.sh', shell=True, cwd=cfd_dir_1, env=cfd_env_1)

    # create new CFD_2 folder
    shutil.copytree('setup_files_GC/fluid', cfd_dir_2)
    for j in range(len(mesh_list)):
        if j != i: remove(cfd_dir_2 + "/" + mesh_list[j] + '_fluid.msh')
    with open(join(cfd_dir_2, 'case_2.jou')) as infile:
        with open(join(cfd_dir_2, 'Case_2.jou'), 'w') as outfile:
            for line in infile:
                line = line.replace('|MESH|', mesh)
                outfile.write(line)
    remove(join(cfd_dir_2, 'case_2.jou'))
    cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
    subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)

    # copy run_simulation.py script and parameters.json to mesh specific directories
    shutil.copy('../run_simulation.py', './' + mesh + '/')
    shutil.copy('../parameters.json', './' + mesh + '/')