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

solver_1 = 'cht_fluent.cht_v2023R1'
solver_2 = 'cht_fluent.cht_v2023R1'

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