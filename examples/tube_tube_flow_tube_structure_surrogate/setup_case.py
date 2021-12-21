import shutil
import subprocess

from coconut import tools

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
shutil.copytree('../setup_files/tube/tube_flow', cfd_dir)

# create new CSM folder
shutil.copytree('../setup_files/tube/tube_structure', csm_dir)

# create new CFD surrogate folder
shutil.copytree('../setup_files/tube/tube_flow', cfd_sur_dir)

# create new CSM surrogate folder
shutil.copytree('../setup_files/tube/tube_structure', csm_sur_dir)
