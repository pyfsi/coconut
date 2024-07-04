import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

round_two = True
pre = './qn_test/'

solver_1 = 'fluent.v2023R1'
solver_2 = 'fluent.v2023R1'

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
