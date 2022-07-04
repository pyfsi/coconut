import shutil
import subprocess

from coconut import tools

cfd_solver = 'fluent.v2019R3'
csm_solver = 'abaqus.v614'

cfd_dir = './CFD'
csm_dir = './CSM'
cfd_sur_dir = './CFD_surrogate'
csm_sur_dir = './CSM_surrogate'

# copy run_simulation.py script to main directory
shutil.copy('../setup_files/run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)
shutil.rmtree(cfd_sur_dir, ignore_errors=True)
shutil.rmtree(csm_sur_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/tube/fluent2d', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call('./setup_fluent2d.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree('../setup_files/tube/abaqus2d', csm_dir)
csm_env = tools.get_solver_env(csm_solver, csm_dir)
subprocess.check_call('./setup_abaqus2d.sh', shell=True, cwd=csm_dir, env=csm_env)

# create new CFD surrogate folder
shutil.copytree('../setup_files/tube/tube_flow', cfd_sur_dir)

# create new CSM surrogate folder
shutil.copytree('../setup_files/tube/tube_structure', csm_sur_dir)

