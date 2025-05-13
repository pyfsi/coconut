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

try:
    subprocess.check_call(['pkill', '-9', 'apip-standalone'])
except subprocess.CalledProcessError as e:
    if e.returncode == 1:
        print("Warning: No 'apip-standalone' process found to kill. Continuing cleanup.")
    else:
        print(f"Error killing 'apip-standalone': {e}")
        # Potentially raise the exception again if this error is critical
        # raise