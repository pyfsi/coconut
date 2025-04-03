import glob
import shutil
import subprocess
from os.path import join

from coconut import tools

cfd_dir_1 = 'CFD_1'
cfd_dir_2 = 'CFD_2'

# clean up Fluent
if glob.glob(join(cfd_dir_1, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_1)
if glob.glob(join(cfd_dir_2, 'cleanup-fluent*')):
    subprocess.check_call(f'sh cleanup-fluent*', shell=True, cwd=cfd_dir_2)

subprocess.check_call(f'kill -9 apip-standalone', shell=True)

# clean working directories
shutil.rmtree(cfd_dir_1, ignore_errors=True)
shutil.rmtree(cfd_dir_2, ignore_errors=True)
