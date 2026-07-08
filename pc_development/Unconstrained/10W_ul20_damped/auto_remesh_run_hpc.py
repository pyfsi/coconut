import os
import shutil
import subprocess
import time
import re

# === USER SETTINGS ===
# Local Cluster Config
sim_setup_machine = 'cfdclu13'  # Machine for setup/data processing
sim_process_machine = 'cfdclu39'  # Machine for data processing
remesh_machine = 'cfdclu13'  # Machine for remeshing

# HPC Config (UGent)
hpc_user = 'vsc43960'
hpc_host = 'login.hpc.ugent.be'
hpc_cluster_module = 'cluster/shinx'  # The module swap required

# Job Scripts
hpc_script_name = 'job_shinx.pbs'  # Script for the FIRST run
hpc_script_restart = 'job_shinx_restart.pbs'  # Script for subsequent runs

# Paths
local_base_dir = os.getcwd()
hpc_base_dir = '/scratch/gent/439/vsc43960/coconut/pc_development/Unconstrained/10W_ul20_damped'

# Workflow Config
num_runs = 550                   # Total number of runs
auto_restart = 0                # Restart with the remeshing script
manual_restart = 19               # New run started with manual remesh

# Instructions in case of restart:
# - Auto restart: remove general pickle file and the node file for the solid
# - Manual restart:
#   - Use the solid node file in DM_converter and move the converted file to a dedicated folder in the setup_file folder
#   - Use the converted node file in Workbench to create the new mesh, use the existing 'remesh' workflow, export the mesh to the dedicated folder created in step 1
#   - Remove the remesh folder in CFD_2 and give the location of the remesh folder in setup_files below in "manual_mesh_local"
#   - Remove the node file in CFD_1 after use
#   - Remove all simulation files HPC-side

folders_to_clean = ["./CFD_1", "./CFD_2"]
files_to_remove = ["fluent.log", "log", "report-file.out", "rigid-body-report-file.out"]

# Path in case of manual remesh
manual_mesh_local = "./setup_files/remesh_19/remesh.msh"
manual_mesh = os.path.join(local_base_dir, manual_mesh_local)

# Sync Configuration
folders_to_sync = ["CFD_1", "CFD_2"]
# These files are Critical: they determine if the run was actually successful
files_to_sync_down = ["case_results.pickle", "case_restart_ts100.pickle"]

# ==========================================
#      HPC INTERACTION FUNCTIONS
# ==========================================

def run_ssh_command(command, host, user, use_hpc_env=False):
    """Executes a command via SSH with a 30s Timeout to prevent hanging."""
    if use_hpc_env:
        full_cmd = (
            f"ssh {user}@{host} "
            f"\"bash -l -c 'module swap {hpc_cluster_module}; cd {hpc_base_dir}; {command}'\""
        )
    else:
        full_cmd = f"ssh -tt {host} \"bash -l -c 'cd {local_base_dir}; {command}'\""

    try:
        # ADDED TIMEOUT: If SSH takes > 30s, kill it and retry later.
        result = subprocess.run(
            full_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            check=False,
            timeout=30
        )

        if result.returncode != 0:
            return f"ERROR: {result.stderr}"

        return result.stdout.strip()

    except subprocess.TimeoutExpired:
        print(f"\n[Warning] SSH command timed out after 30s.")
        return "ERROR: Timeout"


def countdown(t):
    """Prints a live countdown timer."""
    while t > 0:
        mins, secs = divmod(t, 60)
        timer = '{:02d}:{:02d}'.format(mins, secs)
        print(f"\r[WAIT] Next check in: {timer}   ", end="", flush=True)
        time.sleep(1)
        t -= 1
    print("\r[CHECK] Checking status now...      ", end="", flush=True)


def check_all_files_exist_on_hpc(file_list):
    """Checks if ALL files in the list exist on the HPC."""
    checks = " && ".join([f"[ -f {f} ]" for f in file_list])
    cmd = f"{checks} && echo 'ALL_FOUND' || echo 'MISSING'"
    output = run_ssh_command(cmd, hpc_host, hpc_user, use_hpc_env=True)
    return "ALL_FOUND" in output


def delete_files_on_hpc(file_list):
    """Deletes specific files on the HPC."""
    print(f"[HPC] Cleaning up old status files: {file_list}")
    files_str = " ".join(file_list)
    cmd = f"rm -f {files_str}"
    run_ssh_command(cmd, hpc_host, hpc_user, use_hpc_env=True)


def get_job_status(job_id):
    """Checks job status. Returns 'R', 'Q', 'C', 'E' or None (if missing)."""
    clean_id = job_id.split('.')[0]
    # grep returns exit code 1 if not found, so we add || true
    cmd = f"qstat -a | grep {clean_id} || true"
    output = run_ssh_command(cmd, hpc_host, hpc_user, use_hpc_env=True)

    if "ERROR" in output:
        return "SSH_ERROR"

    if not output.strip():
        return None  # Missing from queue

    parts = output.split()
    if len(parts) >= 6:
        return parts[-2]  # Return status column
    return "UNKNOWN"


def monitor_job_completion(job_id, files_to_check):
    """
    Intelligently monitors job with a visual countdown.
    Keeps a history of status checks.
    """
    print(f"\n[HPC] Monitoring Job {job_id}...")

    missing_strikes = 0
    max_strikes = 6  # If job is gone AND files missing for 30 mins -> Error

    while True:
        # 1. Check Job Status
        status = get_job_status(job_id)
        current_time = time.strftime("%H:%M:%S")

        # 2. Handle Active States
        if status in ['R', 'Q', 'H']:
            # We print with a newline so this message PERSISTS in your log
            print(f"\r[{current_time}] Job Status: {status}. Waiting 5 min...")
            missing_strikes = 0
            countdown(300)
            continue

        if status == "SSH_ERROR":
            print(f"\r[{current_time}] SSH Connection glitch. Retrying in 60s...")
            countdown(60)
            continue

        # 3. Handle Completion States (Missing, Completed, Exiting)
        # Clear the countdown line
        print(f"\r[{current_time}] Job status is '{status}' (or missing). Checking files...   ")

        if check_all_files_exist_on_hpc(files_to_check):
            print(f"  -> SUCCESS: Critical files found on HPC.")
            return True

            # 4. Job is gone, but files are MISSING
        missing_strikes += 1
        print(f"  -> Job is finished/gone, but files are missing! (Strike {missing_strikes}/{max_strikes})")

        if missing_strikes >= max_strikes:
            raise FileNotFoundError("Job finished/crashed and critical files never appeared on HPC.")

        print("  -> Waiting 5 min for filesystem sync or job recovery...")
        countdown(300)


def rsync_to_hpc():
    """Uploads entire folders from Local -> HPC."""
    print("\n[SYNC UP] syncing Local -> HPC...")
    for folder in folders_to_sync:
        local_path = os.path.join(local_base_dir, folder)
        if os.path.exists(local_path):
            cmd = f"rsync -avz --delete {local_path} {hpc_user}@{hpc_host}:{hpc_base_dir}/"
            print(f"  -> Uploading {folder}...")
            subprocess.run(cmd, shell=True, check=True)
        else:
            print(f"  [WARNING] Local folder {folder} not found, skipping upload.")


def rsync_from_hpc():
    """Downloads entire folders AND specific files from HPC -> Local."""
    print("\n[SYNC DOWN] syncing HPC -> Local...")

    # 1. Sync Folders
    for folder in folders_to_sync:
        remote_path = f"{hpc_base_dir}/{folder}"
        cmd = f"rsync -avz --delete {hpc_user}@{hpc_host}:{remote_path} {local_base_dir}/"
        print(f"  -> Downloading folder {folder}...")
        subprocess.run(cmd, shell=True, check=True)

    # 2. Sync Specific Files (Pickles)
    print("  -> Downloading specific restart files...")
    file_sources = [f"{hpc_user}@{hpc_host}:{hpc_base_dir}/{f}" for f in files_to_sync_down]
    sources_str = " ".join(file_sources)

    cmd = f"rsync -avz {sources_str} {local_base_dir}/"
    subprocess.run(cmd, shell=True, check=True)


def submit_hpc_job(script_name):
    """Submits the PBS job and returns the Job ID."""
    print(f"[HPC] Submitting {script_name}...")
    output = run_ssh_command(f"qsub {script_name}", hpc_host, hpc_user, use_hpc_env=True)

    match = re.search(r"(\d+)", output)
    if match:
        job_id = match.group(1)
        print(f"[HPC] Job submitted successfully. ID: {job_id}")
        return job_id
    else:
        raise ValueError(f"Could not parse Job ID from output: {output}")


# ==========================================
#      LOCAL HELPERS
# ==========================================

def run_remote_local(script, cluster):
    """Run a Python script on the local cluster (cfdclu machines)."""
    cmd = (
        f"ssh -tt {cluster} "
        f"\"bash -l -c 'module load Anaconda3-python && cd {local_base_dir} && python3 -u {script}'\""
    )
    print(f"[LOCAL:{cluster}] Running: {script}")
    subprocess.run(cmd, shell=True, check=True)


def clean_folder(path):
    """Clean up specified files and directories inside a folder."""
    for root, dirs, files in os.walk(path):
        for f in files:
            if f in files_to_remove or f.startswith("case_timestep"):
                try:
                    os.remove(os.path.join(root, f))
                except OSError:
                    pass
        if "create_mesh" in dirs:
            shutil.rmtree(os.path.join(root, "create_mesh"))


# ==========================================
#      MAIN WORKFLOW
# ==========================================

def main():
    print("=== Starting Hybrid HPC FSI Workflow ===")
    restart = max([manual_restart, auto_restart])

    # Step 1: Initial Setup (Local)
    if not restart:
        run_remote_local("setup_case.py", sim_setup_machine)
        run_remote_local("solid_interface.py", sim_process_machine)

    # === MAIN LOOP ===
    for i in range(max(0, restart - 1), num_runs):
        print(f"\n{'=' * 30}")
        print(f"   STARTING RUN {i + 1}/{num_runs}")
        print(f"{'=' * 30}")

        if i == manual_restart - 1:
            run_remote_local(f"remesh_man.py {manual_mesh}", remesh_machine)
            run_remote_local("solid_interface.py", sim_process_machine)
        if i == restart - 1:
            for folder in folders_to_clean:
                clean_folder(folder)

        # ----------------------------------------
        # PHASE A: PREPARE & SYNC (Local -> HPC)
        # ----------------------------------------
        rsync_to_hpc()

        # ----------------------------------------
        # PHASE B: SIMULATION (HPC)
        # ----------------------------------------
        # Select correct script based on run number
        current_script = hpc_script_name if i == 0 else hpc_script_restart

        job_id = submit_hpc_job(current_script)

        # Monitor until job is done AND files exist
        monitor_job_completion(job_id, files_to_sync_down)

        # ----------------------------------------
        # PHASE C: RETRIEVE & CLEANUP (HPC -> Local)
        # ----------------------------------------
        # 1. Download results
        rsync_from_hpc()

        # 2. DELETE the critical status files from HPC
        # This prevents the NEXT run from seeing old files and thinking it's done instantly.
        delete_files_on_hpc(files_to_sync_down)

        # ----------------------------------------
        # PHASE D: REMESHING (Local)
        # ----------------------------------------
        print(f"\n[LOCAL] Processing data and Remeshing for next step...")

        # 1. Collect Data (Extract forces/displacements)
        run_remote_local("collect_data.py", sim_process_machine)

        # 2. Remesh (Create new geometry/mesh for NEXT run)
        if i < num_runs - 1:
            run_remote_local("remesh.py", remesh_machine)

            # 3. Interface Translation
            run_remote_local("solid_interface.py", sim_process_machine)

            # 4. Cleanup to save space
            for folder in folders_to_clean:
                clean_folder(folder)

    print("\n=== Workflow Completed Successfully ===")


if __name__ == "__main__":
    main()