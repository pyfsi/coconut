import shutil
import subprocess

from coconut import tools

cfd_dir = './CFD'
csm_dir = './CSM'

# copy run_simulation.py script to main directory
shutil.copy('../setup_files/run_simulation.py', './')

# clean working directories
shutil.rmtree(cfd_dir, ignore_errors=True)
shutil.rmtree(csm_dir, ignore_errors=True)

# create new CFD folder
shutil.copytree('../setup_files/tube_flow', cfd_dir)

# create new CSM folder
shutil.copytree('../setup_files/tube_structure', csm_dir)