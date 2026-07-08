import glob
import shutil
import subprocess
import os
from os.path import join

from coconut import tools

solver_2 = 'pc_fluent_liquid_rb.v2024R2'
cfd_dir_1 = 'CFD_1'
cfd_dir_2 = 'CFD_2'

# clean up Fluent
if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

# clean working directories
shutil.rmtree(cfd_dir_1, ignore_errors=True)
shutil.rmtree(cfd_dir_2, ignore_errors=True)

# create new solid folder
shutil.copytree('setup_files/solid', cfd_dir_1)

# Create empty liquid folder
os.makedirs(cfd_dir_2, exist_ok=True)

# Define items to copy (all treated as files)
items_to_copy = ['case_2_flow.jou', 'eicosane', 'setup_fluent_2.sh', 'background_v4.msh', 'component_v10.msh']
src_liquid = join('setup_files', 'liquid')

for item in items_to_copy:
    src_path = join(src_liquid, item)
    dst_path = join(cfd_dir_2, item)

    if os.path.exists(src_path):
        shutil.copy2(src_path, dst_path)
    else:
        print(f"Warning: {item} not found in {src_liquid}")

cfd_env_2 = tools.get_solver_env(solver_2, cfd_dir_2)
subprocess.check_call('./setup_fluent_2.sh', shell=True, cwd=cfd_dir_2, env=cfd_env_2)