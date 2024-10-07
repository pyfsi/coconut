import glob
import shutil
import subprocess
from os.path import join

from coconut import tools

solver_1 = 'cht_fluent.v2023R1'
solver_2 = 'cht_fluent.v2023R1'
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

# create new CFD folder
shutil.copytree('setup_files/solid', cfd_dir_1)
cfd_env_1 = tools.get_solver_env(solver_1, cfd_dir_1)
subprocess.check_call('./setup_fluent_1.sh', shell=True, cwd=cfd_dir_1, env=cfd_env_1)

# create new CSM folder
shutil.copytree('setup_files/liquid', cfd_dir_2)
cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)
