import os
import shutil
import subprocess

# === USER SETTINGS ===
simulation_machine = 'cfdclu33'  # Machine for setup and simulations
remesh_machine = 'cfdclu13'      # Machine for remeshing
num_runs = 9                     # Total number of runs
run_of_tri = 6                   # Run number for switch to tri cells in solid
run_of_smooth = 7                # Run number for switch to non-remeshing in liquid

# Paths to clean
folders_to_clean = ["./CFD_1", "./CFD_2"]

# Files and directories to remove
files_to_remove = ["fluent.log", "log"]
case_pattern = "case_timestep"
dir_to_remove = "create_mesh"

# Get current working directory dynamically
base_dir = os.getcwd()


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
    # Step 1: Setup (on cfdclu33)
    # run_remote("setup_case.py", simulation_machine)

    # Step 2: Initial simulation (on cfdclu33)
    print(f"\n=== Run 1/{num_runs} ===")
    # run_remote("run_simulation.py", simulation_machine) # 14000 time steps

    # Loop over runs (except last one) --> CHANGE BACK, REMOVE THE 1 IN RANGE!!!
    for i in range(3, num_runs - 1):
        # Step 3: Remesh (on cfdclu13)
        run_remote(f"remesh.py {run_of_tri} {run_of_smooth}", remesh_machine) # PASS HERE THE TWO PARAMETERS

        # Step 4: Cleanup (local)
        for folder in folders_to_clean:
            clean_folder(folder, remove_create_mesh=(i == 0))

        # Step 5: Restart simulation (on cfdclu33)
        print(f"\n=== Run {i + 2}/{num_runs} ===")
        run_remote("run_simulation_restart.py", simulation_machine) # 8 x 7200 time steps

    # After loop
    print("\n=== Final Remesh and Cleanup ===")
    run_remote(f"remesh.py {run_of_tri} {run_of_smooth}", remesh_machine)
    for folder in folders_to_clean:
        clean_folder(folder, remove_create_mesh=False)


if __name__ == "__main__":
    main()
