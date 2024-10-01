import glob
import shutil
import subprocess
from os.path import join

from coconut import tools

solver_1 = 'pc_fluent.pc_v2023R1'
solver_2 = 'pc_fluent.pc_v2023R1'
cfd_dir_1 = './CFD_1'
cfd_dir_2 = './CFD_2'

# copy run_simulation.py script to main directory (NOT NECESSARY HERE)
# shutil.copy('../../run_simulation.py', './')

# clean up Fluent
if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

# clean working directories
shutil.rmtree(cfd_dir_1, ignore_errors=True)
shutil.rmtree(cfd_dir_2, ignore_errors=True)
