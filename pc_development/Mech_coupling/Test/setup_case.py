import glob
import shutil
import subprocess
from os.path import join

from coconut import tools

solver_2 = 'pc_fluent.v2024R2'
cfd_dir_2 = './CFD_2'

# copy run_simulation.py script to main directory (NOT NECESSARY HERE)
# shutil.copy('../../run_simulation.py', './')

# clean up Fluent
if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

# clean working directories
shutil.rmtree(cfd_dir_2, ignore_errors=True)

# create new liquid phase CFD folder
shutil.copytree('setup_files/liquid', cfd_dir_2)
cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)
