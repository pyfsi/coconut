import shutil
import subprocess
import os
from os.path import join

from coconut import tools

cfd_solver = 'fluent.v2019R3'
csm_solver = 'abaqus.v614'
cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../../run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/fluent3d', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call('./setup_fluent3d.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree('../setup_files/abaqus3d', csm_dir)
csm_env = tools.get_solver_env(csm_solver, csm_dir)
subprocess.check_call('./setup_abaqus3d.sh', shell=True, cwd=csm_dir, env=csm_env)
