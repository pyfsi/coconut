import shutil
import subprocess

from coconut import tools

cfd_solver = 'fluent.v2021R2'
csm_solver = 'kratos_structure.v91'
cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../../run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/fluent2d', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call('./setup_fluent2d.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree('../setup_files/kratos_structure2d', csm_dir)
