import glob
import shutil
import subprocess
import os
from os.path import join
from os import remove

from coconut import tools

debug = True
reverse = True

if reverse:
    pre = './reverse/'
else:
    pre = './'

solver_1 = 'fluent.v2023R1'
solver_2 = 'fluent.v2023R1'

if debug:
    solver_order = ['TFFB']
    methods = ['relaxation', 'iqni']
else:
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
        print(pre + order + '_' + method + ' removed!')

if not debug:
    if reverse:
        shutil.rmtree(pre, ignore_errors=True)