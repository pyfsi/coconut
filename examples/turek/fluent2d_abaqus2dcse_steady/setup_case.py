import glob
import shutil
import subprocess
from os.path import join

from coconut import tools

cfd_solver = 'fluent.v2024R2'
csm_solver = 'abaqus_cse.v2024'
cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../../run_simulation.py', './')

# clean up Fluent
if glob.glob(join(cfd_dir, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir)

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/fluent2d_steady', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call('./setup_fluent2d_steady.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree('../setup_files/abaqus2d_steadycse', csm_dir)
csm_env = tools.get_solver_env(csm_solver, csm_dir)
