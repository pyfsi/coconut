import subprocess
import yaml
import json
import os
import time
import re


def load_config():
    with open('config.yaml', 'r') as file:
        return yaml.safe_load(file)


def load_state():
    if os.path.exists('state.json'):
        with open('state.json', 'r') as file:
            return json.load(file)
    return {"current_run": 0}


def save_state(state):
    with open('state.json', 'w') as file:
        json.dump(state, file)


# ==========================================
#      LOCAL CLUSTER (cfdclu) FUNCTIONS
# ==========================================

def execute_local_node(script_name, cluster_name, base_dir, extra_args=""):
    """Executes a Python script on the local cluster (e.g., cfdclu13)."""
    command = f"ssh -tt {cluster_name} \"bash -l -c 'module load Anaconda3-python && cd {base_dir} && python3 -u src/{script_name} {extra_args}'\""
    print(f"[LOCAL:{cluster_name}] Running: {script_name} {extra_args}")
    subprocess.run(command, shell=True, check=True)


def clean_local_fluent(cluster_name, base_dir, cfd2_dir):
    """Kills orphan Fluent and API processes on local clusters."""
    print(f"\n[CLEANUP] Executing local Fluent cleanup on {cluster_name}...")
    command = (
        f"ssh -tt {cluster_name} \"bash -l -c '"
        f"cd {base_dir}/{cfd2_dir} && "
        f"if ls cleanup-fluent* 1> /dev/null 2>&1; then sh cleanup-fluent*; fi; "
        f"pkill -9 apip-standalone'\""
    )
    try:
        subprocess.run(command, shell=True, check=False)
        print("[CLEANUP] Local cleanup finished successfully.\n")
    except Exception as e:
        print(f"[CLEANUP] Error during local cleanup: {e}\n")


def prepare_parameters(config, is_restart):
    """Dynamically updates parameters.json for initial or restart runs."""
    param_file = config['coconut']['base_parameter_file']

    # Ensure the file exists
    if not os.path.exists(param_file):
        raise FileNotFoundError(f"Missing base parameter file: {param_file}")

    with open(param_file, 'r') as f:
        params = json.load(f)

    # Grab the target settings from the config based on the run phase
    if is_restart:
        mode_settings = config['coconut']['restart_run']
        mode_name = "RESTART"
    else:
        mode_settings = config['coconut']['initial_run']
        mode_name = "INITIAL"

    # Find the liquid wrapper and inject the new settings
    for wrapper in params['coupled_solver']['solver_wrappers']:
        if 'pc_fluent_liquid' in wrapper['type']:
            wrapper['settings']['case_file'] = mode_settings['case_file']
            wrapper['settings']['RB']['restart_rb'] = mode_settings['restart_rb']
            break

    # Overwrite the file on disk
    with open(param_file, 'w') as f:
        json.dump(params, f, indent=2)

    print(f"[SETUP] Dynamically configured {param_file} for {mode_name} mode.")

# ==========================================
#      HPC (UGent) FUNCTIONS
# ==========================================

def run_hpc_command(command, config):
    """Executes a command on the HPC via SSH with a 30s timeout."""
    hpc = config['hpc']
    full_cmd = (
        f"ssh {hpc['user']}@{hpc['host']} "
        f"\"bash -l -c 'module swap {hpc['module']}; cd {hpc['base_dir']}; {command}'\""
    )
    try:
        result = subprocess.run(full_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=False, timeout=30)
        if result.returncode != 0:
            return f"ERROR: {result.stderr}"
        return result.stdout.strip()
    except subprocess.TimeoutExpired:
        print(f"\n[Warning] SSH command timed out after 30s.")
        return "ERROR: Timeout"


def rsync_to_hpc(config):
    print("\n[SYNC UP] Syncing Local -> HPC...")
    hpc = config['hpc']
    for folder in config['sync']['folders']:
        local_path = os.path.join(config['directories']['base_dir'], folder)
        if os.path.exists(local_path):
            cmd = f"rsync -avz --delete {local_path} {hpc['user']}@{hpc['host']}:{hpc['base_dir']}/"
            print(f"  -> Uploading {folder}...")
            subprocess.run(cmd, shell=True, check=True)


def rsync_from_hpc(config):
    print("\n[SYNC DOWN] Syncing HPC -> Local...")
    hpc = config['hpc']
    local_base = config['directories']['base_dir']

    for folder in config['sync']['folders']:
        remote_path = f"{hpc['base_dir']}/{folder}"
        cmd = f"rsync -avz --delete {hpc['user']}@{hpc['host']}:{remote_path} {local_base}/"
        print(f"  -> Downloading folder {folder}...")
        subprocess.run(cmd, shell=True, check=True)

    print("  -> Downloading specific critical files...")
    file_sources = [f"{hpc['user']}@{hpc['host']}:{hpc['base_dir']}/{f}" for f in config['sync']['critical_files']]
    cmd = f"rsync -avz {' '.join(file_sources)} {local_base}/"
    subprocess.run(cmd, shell=True, check=True)


def submit_hpc_job(script_name, config):
    print(f"[HPC] Submitting {script_name}...")
    output = run_hpc_command(f"qsub {script_name}", config)
    match = re.search(r"(\d+)", output)
    if match:
        job_id = match.group(1)
        print(f"[HPC] Job submitted successfully. ID: {job_id}")
        return job_id
    raise ValueError(f"Could not parse Job ID from output: {output}")


def cancel_hpc_job(job_id, config):
    if job_id:
        print(f"\n[HPC CLEANUP] Cancelling remote job {job_id}...")
        run_hpc_command(f"qdel {job_id}", config)


def delete_files_on_hpc(config):
    file_list = config['sync']['critical_files']
    print(f"[HPC] Cleaning up old status files: {file_list}")
    run_hpc_command(f"rm -f {' '.join(file_list)}", config)


def countdown(t):
    while t > 0:
        mins, secs = divmod(t, 60)
        print(f"\r[WAIT] Next check in: {mins:02d}:{secs:02d}   ", end="", flush=True)
        time.sleep(1)
        t -= 1
    print("\r[CHECK] Checking status now...      ", end="", flush=True)


def monitor_job_completion(job_id, config):
    print(f"\n[HPC] Monitoring Job {job_id}...")
    missing_strikes = 0
    max_strikes = 6

    while True:
        clean_id = job_id.split('.')[0]
        output = run_hpc_command(f"qstat -a | grep {clean_id} || true", config)
        current_time = time.strftime("%H:%M:%S")

        if "ERROR" in output:
            print(f"\r[{current_time}] SSH Connection glitch. Retrying in 60s...")
            countdown(60)
            continue

        status = output.split()[-2] if output.strip() and len(output.split()) >= 6 else None

        if status in ['R', 'Q', 'H']:
            print(f"\r[{current_time}] Job Status: {status}. Waiting 5 min...\n")
            missing_strikes = 0
            countdown(300)
            continue

        print(f"\r[{current_time}] Job status is '{status}' (or missing). Checking files...   ")
        checks = " && ".join([f"[ -f {f} ]" for f in config['sync']['critical_files']])
        file_check = run_hpc_command(f"{checks} && echo 'ALL_FOUND' || echo 'MISSING'", config)

        if "ALL_FOUND" in file_check:
            print("  -> SUCCESS: Critical files found on HPC.")
            return True

        missing_strikes += 1
        print(f"  -> Job is finished/gone, but files are missing! (Strike {missing_strikes}/{max_strikes})")

        if missing_strikes >= max_strikes:
            raise RuntimeError("Job finished/crashed and critical files never appeared on HPC.")

        print("  -> Waiting 5 min for filesystem sync or job recovery...")
        countdown(300)