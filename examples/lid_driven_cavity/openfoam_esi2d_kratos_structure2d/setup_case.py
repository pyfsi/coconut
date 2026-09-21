import shutil
import subprocess

from coconut import tools

# cfd_solver = 'openfoam_esi.v2312'
cfd_name = 'openfoam_esi'
cfd_solver = f'{cfd_name}.v2312'
csm_name = 'kratos_structure'
csm_solver = f'{csm_name}.v94'
cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../../run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree(f'../setup_files/{cfd_name}2d', cfd_dir)
cfd_env = tools.get_solver_env(cfd_solver, cfd_dir)
subprocess.check_call(f'./setup_{cfd_name}2d.sh', shell=True, cwd=cfd_dir, env=cfd_env)

# create new CSM folder
shutil.copytree(f'../setup_files/{csm_name}2d', csm_dir)
