import shutil
import subprocess

from coconut import tools

cfd_solver = 'fluent.v2019R1'
csm_solver = 'abaqus.v614'
cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../setup_files/run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/turek/fluent2d_steady', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call('./setup_fluent2d_steady.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree('../setup_files/turek/abaqus2d_steady', csm_dir)
csm_env = tools.get_solver_env(csm_solver, csm_dir)
subprocess.check_call('./setup_abaqus2d_steady.sh', shell=True, cwd=csm_dir, env=csm_env)
