import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

round_two = True
reverse = False
fast = False

pre = './qn_test/'

solver_1 = 'cht_fluent.v2023R1'
solver_2 = 'cht_fluent.v2023R1'

# settings
mesh = 'M3'
solid_cores = 1
fluid_cores = 1
conduction = 0.2876 # W/mK
method = 'iqni'
omega = 0.2

if round_two:
    solver_order = ['TFFB', 'FFTB']
    min_sign = [1e-4, 1e-2]
    q_list = [3]
else:
    solver_order = ['TFFB', 'FFTB']
    min_sign = [1e-13, 1e-6]
    q_list = [0, 1, 2, 3, 4]

for i, order in enumerate(solver_order):
    for q in q_list:
        for ms in min_sign:
            cfd_dir_1 = pre + order + '/q' + str(q) + '/' + str(ms) + '/CFD_1'
            cfd_dir_2 = pre + order + '/q' + str(q) + '/' + str(ms) + '/CFD_2'

            # clean up Fluent
            if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
                subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
            if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
                subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

            # clean working directories
            shutil.rmtree(pre + order + '/q' + str(q) + '/' + str(ms), ignore_errors=True)

if not round_two:
    shutil.rmtree(pre, ignore_errors=True)

try:
    os.mkdir(pre)
except:
    print(pre + ' already exists, this is round 2.')

for i, order in enumerate(solver_order):
    try:
        os.mkdir(pre + order)
    except:
        print(pre + order + ' already exists.')
    for q in q_list:
        try:
            os.mkdir(pre + order + '/q' + str(q))
        except:
            print(pre + order + '/q' + str(q) + ' already exists.')
        for ms in min_sign:
            os.mkdir(pre + order + '/q' + str(q) + '/' + str(ms))

            cfd_dir_1 = pre + order + '/q' + str(q) + '/' + str(ms) + '/CFD_1'
            cfd_dir_2 = pre + order + '/q' + str(q) + '/' + str(ms) + '/CFD_2'

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
            shutil.copy('setup_files/run_simulation.py', pre + order + '/q' + str(q) + '/' + str(ms) + '/')
            with open('setup_files/parameters_' + order + '.json') as infile:
                with open(pre + order + '/q' + str(q) + '/' + str(ms) + '/parameters.json', 'w') as outfile:
                    for line in infile:
                        line = line.replace('|MESH|', mesh)
                        line = line.replace('|SOLID_CORES|', str(solid_cores))
                        line = line.replace('|FLUID_CORES|', str(fluid_cores))
                        line = line.replace('|METHOD|', method)
                        line = line.replace('|OMEGA|', str(omega))
                        line = line.replace('|OMEGA_MAX|', str(omega))
                        line = line.replace('|MIN_SIGN|', str(ms))
                        line = line.replace('|Q|', str(q))
                        outfile.write(line)