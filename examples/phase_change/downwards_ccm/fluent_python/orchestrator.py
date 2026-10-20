import sys
import subprocess
import argparse
from src.utils import (
    load_config, load_state, save_state, execute_local_node,
    clean_local_fluent, rsync_to_hpc, rsync_from_hpc, submit_hpc_job,
    monitor_job_completion, delete_files_on_hpc, cancel_hpc_job,
    prepare_parameters
)


def main():
    parser = argparse.ArgumentParser(description="Hybrid HPC FSI Workflow Orchestrator")
    parser.add_argument("--manual-mesh", action="store_true",
                        help="Resume the workflow using manual remeshing for the first run.")
    args = parser.parse_args()

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

    # Guard against accidental re-runs if all runs are finished
    if start_run >= total_runs:
        print(f"\n[INFO] Workflow already completed ({start_run}/{total_runs} runs recorded in state.json).")
        return

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
            handle_crash(0, mac_setup, mac_process, config, current_job_id, e)

    # Phase B: Remeshing Loop
    first_manual_override = args.manual_mesh
    skip_remesh = False

    while state['current_run'] < total_runs:
        i = state['current_run']
        print(f"\n{'=' * 30}")
        print(f"   STARTING RUN {i + 1}/{total_runs}")
        print(f"{'=' * 30}")

        try:
            if first_manual_override:
                print(f"\n[OVERRIDE] Executing initial manual remesh for Run {i + 1}...")
                manual_mesh_path = f"./{config['directories']['manual_interventions']}/new_mesh.msh"
                execute_local_node("remesh_pipeline.py", mac_remesh, base_dir,
                                   extra_args=f"--manual-mesh {manual_mesh_path}")
                execute_local_node("solid_interface.py", mac_process, base_dir)
                prepare_parameters(config, is_restart=True)
                first_manual_override = False

            elif not skip_remesh:
                execute_local_node("collect_data.py", mac_process, base_dir)
                execute_local_node("remesh_pipeline.py", mac_remesh, base_dir)
                execute_local_node("solid_interface.py", mac_process, base_dir)
                prepare_parameters(config, is_restart=True)

            skip_remesh = False

            # Simulation Execution
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
            recovered = handle_crash(i, mac_setup, mac_process, config, current_job_id, e)
            current_job_id = None
            if recovered:
                skip_remesh = True

    # Phase C: Finalisation (Collect last run data and archive final state)
    print("\n" + "=" * 50)
    print("   FINALISING WORKFLOW (POST-PROCESSING & ARCHIVE)")
    print("=" * 50)
    execute_local_node("collect_data.py", mac_process, base_dir)
    execute_local_node("archive_run.py", mac_process, base_dir)

    print("\n=== Workflow Completed Successfully ===")


def handle_crash(run_num, mac_setup, mac_process, config, job_id, error):
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
        base_dir = config['directories']['base_dir']
        manual_dir = config['directories']['manual_interventions']
        manual_mesh_path = f"./{manual_dir}/new_mesh.msh"
        dm_curve_file = f"./{manual_dir}/interface_nodes_designmodeler.txt"

        print("\n" + "=" * 65)
        print("              ENTERING MANUAL RECOVERY MODE")
        print("=" * 65)
        print("\n[STEP 1/3] Converting CFD_1/interface_nodes.csv for ANSYS DesignModeler...")
        execute_local_node("export_designmodeler.py", mac_process, base_dir)

        print("\n[STEP 2/3] DesignModeler curve coordinates ready at:")
        print(f"  -> {dm_curve_file}")
        print("\nACTION REQUIRED:")
        print("  1. Open ANSYS DesignModeler.")
        print("  2. In DesignModeler, create a 3D Curve using the coordinates from:")
        print(f"     {dm_curve_file}")
        print("  3. Build your updated mesh and save the final .msh file to:")
        print(f"     {manual_mesh_path}")
        print("=" * 65)

        input("\n[STEP 3/3] Press ENTER once new_mesh.msh is saved in place...")

        # Inject manual mesh and prepare restart state
        execute_local_node("remesh_pipeline.py", config['machines']['remesh'], base_dir,
                           extra_args=f"--manual-mesh {manual_mesh_path}")
        execute_local_node("solid_interface.py", mac_process, base_dir)
        prepare_parameters(config, is_restart=True)

        return True
    else:
        print("Invalid choice. Defaulting to abort.")
        sys.exit(1)


if __name__ == "__main__":
    main()