import glob
import shutil
import subprocess
import os
from os.path import join

debug = False

solver_1 = 'cht_fluent.v2023R1'
solver_2 = 'cht_fluent.v2023R1'

velocity_list = [2.4, 12, 24, 36] # m/s
velocity_names = ['2-4m_s', '12m_s', '24m_s', '36m_s']
solver_order = ['TFFB', 'FFTB']
methods = ['relaxation', 'aitken', 'iqni']

if debug:
    velocity_list = [2.4]
    velocity_names = ['2-4m_s']

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