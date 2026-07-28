import sys
import subprocess
from src.utils import (
    load_config, load_state, save_state, execute_local_node,
    clean_local_fluent, rsync_to_hpc, rsync_from_hpc, submit_hpc_job,
    monitor_job_completion, delete_files_on_hpc, cancel_hpc_job,
    prepare_parameters
)


def main():
    print("=== Starting Hybrid HPC FSI Workflow ===")
    config = load_config()
    state = load_state()

    mac_setup = config['machines']['local_setup']
    mac_process = config['machines']['local_process']
    mac_remesh = config['machines']['remesh']
    base_dir = config['directories']['base_dir']

    use_hpc = config['hpc']['enabled']
    total_runs = config['simulation']['total_runs']
    start_run = state['current_run']

    # Track Job ID so we can cancel it if the user interrupts
    current_job_id = None

    # Phase A: Initialisation (Run 0)
    if start_run == 0:
        print("\n=== Initialising New Workflow ===")
        execute_local_node("setup_case.py", mac_setup, base_dir)
        execute_local_node("solid_interface.py", mac_process, base_dir)

        # Configure parameters for INITIAL run
        prepare_parameters(config, is_restart=False)

        try:
            if use_hpc:
                rsync_to_hpc(config)
                # Now uses the single unified PBS script
                current_job_id = submit_hpc_job(config['hpc']['script'], config)
                monitor_job_completion(current_job_id, config)
                rsync_from_hpc(config)
                delete_files_on_hpc(config)
            else:
                execute_local_node("run_simulation.py", mac_setup, base_dir)

            # Save state only upon total success
            state['current_run'] = 1
            save_state(state)

        except (KeyboardInterrupt, subprocess.CalledProcessError, RuntimeError) as e:
            handle_crash(0, mac_setup, config, current_job_id, e)

    # Phase B: Remeshing Loop
    for i in range(state['current_run'], total_runs):
        print(f"\n{'=' * 30}")
        print(f"   STARTING RUN {i + 1}/{total_runs}")
        print(f"{'=' * 30}")

        try:
            # 1. Local Data Processing & Remeshing
            execute_local_node("collect_data.py", mac_process, base_dir)
            execute_local_node("remesh_pipeline.py", mac_remesh, base_dir)
            execute_local_node("solid_interface.py", mac_process, base_dir)

            # Configure parameters for RESTART run
            prepare_parameters(config, is_restart=True)

            # 2. Simulation Execution
            if use_hpc:
                rsync_to_hpc(config)
                current_job_id = submit_hpc_job(config['hpc']['script'], config)
                monitor_job_completion(current_job_id, config)
                rsync_from_hpc(config)
                delete_files_on_hpc(config)
            else:
                # Use run_simulation.py for both modes as parameters are handled dynamically
                execute_local_node("run_simulation.py", mac_setup, base_dir)

            # 3. Save State
            state['current_run'] = i + 1
            save_state(state)
            current_job_id = None  # Reset for next loop

        except (KeyboardInterrupt, subprocess.CalledProcessError, RuntimeError) as e:
            handle_crash(i, mac_setup, config, current_job_id, e)
            current_job_id = None

    print("\n=== Workflow Completed Successfully ===")


def handle_crash(run_num, mac_setup, config, job_id, error):
    """Handles manual interrupts and remote crashes gracefully."""
    print(f"\n[INTERRUPT DETECTED] Simulation halted at run {run_num + 1}.")
    print(f"Reason: {error}")

    # 1. Clean up lingering processes
    if config['hpc']['enabled'] and job_id:
        cancel_hpc_job(job_id, config)
    else:
        clean_local_fluent(mac_setup, config['directories']['base_dir'], config['directories']['cfd2'])

    # 2. Prompt for Recovery
    print("\nOptions:")
    print("  [1] Abort the workflow completely.")
    print("  [2] Pause and wait for a manual remesh.")

    choice = input("Enter your choice (1 or 2): ")

    if choice == '1':
        print("Aborting workflow. Goodbye.")
        sys.exit(1)
    elif choice == '2':
        print("\nEntering manual recovery mode.")
        print(f"Please build a new mesh and save it to: ./{config['directories']['manual_interventions']}/new_mesh.msh")
        input("Press ENTER when the new mesh is ready to resume...")

        # Trigger the unified remesh pipeline with the manual mesh argument
        manual_mesh_path = f"./{config['directories']['manual_interventions']}/new_mesh.msh"
        execute_local_node("remesh_pipeline.py", config['machines']['remesh'],
                           config['directories']['base_dir'],
                           extra_args=f"--manual-mesh {manual_mesh_path}")

        # We do not advance the state here so the 'for' loop will retry the same run_num
    else:
        print("Invalid choice. Defaulting to abort.")
        sys.exit(1)


if __name__ == "__main__":
    main()