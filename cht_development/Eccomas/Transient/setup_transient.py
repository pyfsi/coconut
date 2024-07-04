import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

reverse = False
fast = False

if reverse:
    if fast:
        pre = './reverse_fast/'
    else:
        pre = './reverse_slow/'
else:
    if fast:
        pre = './fast/'
    else:
        pre = './slow/'

solver_1 = 'fluent.v2023R1'
solver_2 = 'fluent.v2023R1'

# settings
mesh = 'M3'
solid_cores = 1
fluid_cores = 4
conduction = 0.2876 # W/mK
min_sign = 1e-13
q = 3

# relaxation factors
omega_gs = [0.05, 0.2]
omega_aitken = [0.4, 0.4]
omega_iqni = 0.2

def relax(i, method):
    if method == 'iqni':
        return omega_iqni
    elif method == 'relaxation':
        return omega_gs[i]
    elif method == 'aitken':
        return omega_aitken[i]
    else:
        raise ValueError(f'Method {method} not implemented.')

solver_order = ['TFFB', 'FFTB']
methods = ['relaxation', 'aitken', 'iqni']

for i, order in enumerate(solver_order):
    for k, method in enumerate(methods):

        cfd_dir_1 = pre + order + '_' + method + '/CFD_1'
        cfd_dir_2 = pre + order + '_' + method + '/CFD_2'

        # clean up Fluent
        if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
            subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
        if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
            subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

        # clean working directories
        shutil.rmtree(pre + order + '_' + method, ignore_errors=True)

    shutil.rmtree(pre, ignore_errors=True)
    os.mkdir(pre)

for i, order in enumerate(solver_order):
    for k, method in enumerate(methods):

        os.mkdir(pre + order + '_' + method)

        cfd_dir_1 = pre + order + '_' + method + '/CFD_1'
        cfd_dir_2 = pre + order + '_' + method + '/CFD_2'

        # create new CFD_1 folder
        os.mkdir(cfd_dir_1)
        shutil.copy('setup_files/solid' + "/" + mesh + '_solid.msh', cfd_dir_1)
        shutil.copy('setup_files/solid/setup_fluent_1.sh', cfd_dir_1)
        if order == "FFTB":
            if reverse:
                shutil.copy('setup_files/solid/Q_1_rev.prof', cfd_dir_1)
            else:
                shutil.copy('setup_files/solid/Q_1.prof', cfd_dir_1)
        else:
            if reverse:
                shutil.copy('setup_files/solid/T_1_rev.prof', cfd_dir_1)
            else:
                shutil.copy('setup_files/solid/T_1.prof', cfd_dir_1)
        with open('setup_files/solid/case_1.jou') as infile:
            with open(join(cfd_dir_1, 'case_1.jou'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MESH|', mesh)
                    line = line.replace('|ORDER|', "#t" if order == 'FFTB' else "#f")
                    line = line.replace('|REVERSE|', "#t" if reverse else "#f")
                    line = line.replace('|COND|', str(conduction))
                    outfile.write(line)
        cfd_env_1 = tools.get_solver_env(solver_1, cfd_dir_1)
        subprocess.check_call('./setup_fluent_1.sh', shell=True, cwd=cfd_dir_1, env=cfd_env_1)

        # create new CFD_2 folder
        os.mkdir(cfd_dir_2)
        shutil.copy('setup_files/fluid' + "/" + mesh + '_fluid.msh', cfd_dir_2)
        shutil.copy('setup_files/fluid/setup_fluent_2.sh', cfd_dir_2)
        if order == "FFTB":
            if reverse:
                shutil.copy('setup_files/fluid/T_2_rev.prof', cfd_dir_2)
            else:
                shutil.copy('setup_files/fluid/T_2.prof', cfd_dir_2)
        else:
            if reverse:
                shutil.copy('setup_files/fluid/Q_2_rev.prof', cfd_dir_2)
            else:
                shutil.copy('setup_files/fluid/Q_2.prof', cfd_dir_2)
        if fast:
            journal_file = 'case_2_fast.jou'
        else:
            journal_file = 'case_2_slow.jou'
        with open('setup_files/fluid/' + journal_file) as infile:
            with open(join(cfd_dir_2, 'case_2.jou'), 'w') as outfile:
                for line in infile:
                    line = line.replace('|MESH|', mesh)
                    line = line.replace('|ORDER|', "#t" if order == 'FFTB' else "#f")
                    line = line.replace('|REVERSE|', "#t" if reverse else "#f")
                    outfile.write(line)
        cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
        subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)

        # copy run_simulation.py script and parameters.json to mesh specific directories
        shutil.copy('setup_files/run_simulation.py', pre + order + '_' + method + '/')
        with open('setup_files/parameters_' + order + '.json') as infile:
            with open(pre + order + '_' + method + '/parameters.json', 'w') as outfile:
                for line in infile:
                    line = line.replace('|MESH|', mesh)
                    line = line.replace('|SOLID_CORES|', str(solid_cores))
                    line = line.replace('|FLUID_CORES|', str(fluid_cores))
                    line = line.replace('|METHOD|', method)
                    line = line.replace('|OMEGA|', str(relax(i, method)))
                    line = line.replace('|OMEGA_MAX|', str(relax(i, method)))
                    line = line.replace('|MIN_SIGN|', str(min_sign))
                    line = line.replace('|Q|', str(q))
                    outfile.write(line)