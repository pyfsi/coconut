import glob
import shutil
import subprocess
from os.path import join

solver_1 = 'cht_fluent.v2023R1'
solver_2 = 'cht_fluent.v2023R1'
mesh_list = ['M1', 'M2', 'M3', 'M4']

for i, mesh in enumerate(mesh_list):

    cfd_dir_1 = './' + mesh + '/CFD_1'
    cfd_dir_2 = './' + mesh + '/CFD_2'

    # clean up Fluent
    if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
    if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
        subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

    # clean mesh directories
    shutil.rmtree('./' + mesh, ignore_errors=True)