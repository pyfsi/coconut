import os
import shutil
import subprocess

# === USER SETTINGS ===
simulation_machine = 'cfdclu33'  # Machine for setup and simulations
remesh_machine = 'cfdclu13'      # Machine for remeshing
num_runs = 550                   # Total number of runs
auto_restart = 0
manual_restart = 540             # New run started with manual remesh

# Paths to clean
folders_to_clean = ["./CFD_1", "./CFD_2"]

#Peebsieisdeleukstevrouwvandewereld!
# Files and directories to remove
files_to_remove = ["fluent.log", "log", "report-file.out", "rigid-body-report-file.out"]
case_pattern = "case_timestep"
dir_to_remove = "create_mesh"

# Get current working directory dynamically
base_dir = os.getcwd()

# Path in case of manual remesh
manual_mesh_local = "./setup_files/remesh_540/remesh.msh"
manual_mesh = os.path.join(base_dir, manual_mesh_local)

# Perform restart checks
if manual_restart:
    if auto_restart:
        raise ValueError("If manual restart is active, auto restart should be inactive.")
    if not os.path.exists(manual_mesh):
        raise FileNotFoundError(f"Manual restart is active, but the mesh '{manual_mesh}' was not found.")

def run_local(command):
    """Run a command locally."""
    print(f"[LOCAL] Running: {command}")
    subprocess.run(command, shell=True, check=True)


def run_remote(script, cluster):
    """Run a Python script on a remote cluster via SSH with live output."""
    ssh_command = (
        f"ssh -tt {cluster} "
        f"\"bash -l -c 'module load Anaconda3-python && cd {base_dir} && python3 -u {script}'\""
    )
    print(f"[REMOTE:{cluster}] Running: {script}")
    subprocess.run(ssh_command, shell=True, check=True)



def clean_folder(path, remove_create_mesh=False):
    """Clean up specified files and directories inside a folder."""
    for root, dirs, files in os.walk(path):
        # Remove unwanted files
        for f in files:
            if f in files_to_remove or f.startswith(case_pattern):
                file_path = os.path.join(root, f)
                print(f"Removing file: {file_path}")
                os.remove(file_path)

        # Remove create_mesh directory if required
        if remove_create_mesh and dir_to_remove in dirs:
            dir_path = os.path.join(root, dir_to_remove)
            print(f"Removing directory: {dir_path}")
            shutil.rmtree(dir_path)


def main():
    restart = max([manual_restart, auto_restart])

    # Step 1: Setup
    if not restart:
        run_remote("setup_case.py", simulation_machine)
        # Translate interface to solid side
        run_remote("solid_interface.py", simulation_machine)

    # Step 2: Initial simulation
    if restart <= 1:
        print(f"\n=== Run 1/{num_runs} ===")
        run_remote("run_simulation.py", simulation_machine)

    # Loop over runs
    for i in range(max(0, restart - 2), num_runs - 1):
        # Step 3: Remesh
        if i > restart - 2:
            run_remote("collect_data.py", simulation_machine)
        if i == manual_restart - 2:
            run_remote(f"remesh_man.py {manual_mesh}", remesh_machine)
        elif i > restart - 2:
            run_remote(f"remesh.py", remesh_machine)
        # Translate interface to solid side
        run_remote("solid_interface.py", simulation_machine)

        # Step 4: Cleanup
        for folder in folders_to_clean:
            clean_folder(folder, remove_create_mesh=(i == 0))

        # Step 5: Restart simulation
        print(f"\n=== Run {i + 2}/{num_runs} ===")
        run_remote("run_simulation_restart.py", simulation_machine)

    # After loop
    print("\n=== Final Remesh and Cleanup ===")
    run_remote("collect_data.py", simulation_machine)
    run_remote(f"remesh.py", remesh_machine)
    for folder in folders_to_clean:
        clean_folder(folder, remove_create_mesh=False)


if __name__ == "__main__":
    main()
