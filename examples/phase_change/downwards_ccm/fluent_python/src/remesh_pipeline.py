import os
import shutil
import glob
import re
import subprocess
import argparse
import sys


def run_cfd_workflow(manual_mesh_path=None):
    """
    Automates the CFD simulation workflow including file cleanup,
    directory management, and Fluent execution.
    Handles both automated and manual remeshing paths.
    """
    if manual_mesh_path:
        print(f"Starting CFD workflow in MANUAL mode using mesh: {manual_mesh_path}")
    else:
        print("Starting CFD workflow in AUTOMATED mode...")

    # --- Define main directories and files ---
    current_dir = os.getcwd()
    cfd1_dir = os.path.join(current_dir, 'CFD_1')
    cfd2_dir = os.path.join(current_dir, 'CFD_2')
    setup_files_solid_dir = os.path.join(current_dir, 'setup_files', 'solid')
    setup_files_liquid_dir = os.path.join(current_dir, 'setup_files', 'liquid')
    setup_files_solid_remesh_dir = os.path.join(setup_files_solid_dir, 'create_mesh')
    setup_files_liquid_remesh_dir = os.path.join(setup_files_liquid_dir, 'create_mesh')

    # Ensure CFD_1 and CFD_2 directories exist
    if not os.path.isdir(cfd1_dir) or not os.path.isdir(cfd2_dir):
        print(f"Error: 'CFD_1' or 'CFD_2' directory not found in '{current_dir}'. Please ensure they exist.")
        return

    # --- Step 1: Remove specified files in CFD_1 and CFD_2 ---
    print("\nStep 1: Cleaning up CFD_1 and CFD_2 directories...")
    files_to_remove_patterns = [
        '*.coco', '*.trn', 'cleanup-*.sh', 'displacement_timestep*',
        'heat_flux_timestep*', 'nodes_pc*', 'nodes_update*',
        'RB_update*', 'rigid_body_timestep*'
    ]

    for cfd_dir in [cfd1_dir, cfd2_dir]:
        print(f"  Cleaning files in {cfd_dir}...")
        for pattern in files_to_remove_patterns:
            for file_path in glob.glob(os.path.join(cfd_dir, pattern)):
                try:
                    os.remove(file_path)
                    print(f"    Removed: {file_path}")
                except OSError as e:
                    print(f"    Error removing {file_path}: {e}")

    # --- Step 2: Determine run number and evaluate archive paths ---
    print("\nStep 2: Determining run number and evaluating archive paths...")
    existing_run_folders = glob.glob(os.path.join(current_dir, 'run_*'))
    max_run_num = 0
    for folder in existing_run_folders:
        match = re.match(r'run_(\d+)', os.path.basename(folder))
        if match:
            max_run_num = max(max_run_num, int(match.group(1)))

    if manual_mesh_path:
        # In manual recovery mode, do not archive the crashed attempt into a new run folder.
        # Instead, retrieve clean restart data from the last successfully completed run archive.
        if max_run_num > 0:
            search_dir = os.path.join(current_dir, f'run_{max_run_num}', 'CFD_2')
            print(f"  [MANUAL MODE] Restoring from last completed archive: run_{max_run_num}")
        else:
            search_dir = cfd2_dir
            print("  [MANUAL MODE] No previous run archives found, using active CFD_2 directory")
        previous_run_folder_path = None
        previous_run_folder_name = None
        is_first_run = False
    else:
        past_run_num = max_run_num + 1
        previous_run_folder_name = f'run_{past_run_num}'
        previous_run_folder_path = os.path.join(current_dir, previous_run_folder_name)
        is_first_run = (past_run_num == 1)
        search_dir = cfd2_dir

        try:
            os.makedirs(previous_run_folder_path, exist_ok=True)
            print(f"  Created run folder: {previous_run_folder_path}")
        except OSError as e:
            print(f"  Error creating folder {previous_run_folder_path}: {e}")
            return

    # --- Step 3: Copy CFD folders and results to the archive folder (Automated Mode only) ---
    if not manual_mesh_path:
        print(f"\nStep 3: Copying CFD folders and results to {previous_run_folder_name}...")
        files_to_copy_to_run = ['case_results.pickle', 'case_restart*', 'restart_rb_*.pickle']

        for src_cfd_dir in [cfd1_dir, cfd2_dir]:
            dest_cfd_dir = os.path.join(previous_run_folder_path, os.path.basename(src_cfd_dir))
            try:
                if os.path.exists(dest_cfd_dir):
                    shutil.rmtree(dest_cfd_dir)
                shutil.copytree(src_cfd_dir, dest_cfd_dir)
                print(f"  Copied '{os.path.basename(src_cfd_dir)}' to '{previous_run_folder_name}'")
            except shutil.Error as e:
                print(f"  Error copying directory {src_cfd_dir} to {dest_cfd_dir}: {e}")

        for pattern in files_to_copy_to_run:
            for file_path in glob.glob(os.path.join(current_dir, pattern)):
                try:
                    shutil.copy2(file_path, previous_run_folder_path)
                    print(f"  Copied '{os.path.basename(file_path)}' to '{previous_run_folder_name}'")
                except Exception as e:
                    print(f"  Error copying {file_path}: {e}")
    else:
        print("\nStep 3: Skipping archive creation (MANUAL MODE: preserving previous successful archive)...")

    # --- Step 4: Remove local results files (keeping restart_rb_*.pickle) ---
    print("\nStep 4: Removing local results files...")
    files_to_remove_local = ['case_results.pickle', 'case_restart*']
    for pattern in files_to_remove_local:
        for file_path in glob.glob(os.path.join(current_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"  Removed: {file_path}")
            except OSError as e:
                print(f"  Error removing {file_path}: {e}")

    # Remove report files from CFD_2 so they do not append indefinitely
    files_to_remove_cfd2_reports = ['*report*.out']
    for pattern in files_to_remove_cfd2_reports:
        for file_path in glob.glob(os.path.join(cfd2_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"  Removed: {file_path}")
            except OSError:
                pass

    # --- Step 5 & 6: Copy 'step_1' and 'step_2' jou files ---
    print("\nStep 5 & 6: Copying Fluent journal files...")
    try:
        shutil.copy2(os.path.join(setup_files_liquid_dir, 'step_1.jou'), cfd2_dir)
        print(f"  Copied 'step_1.jou' to '{cfd2_dir}'")
        shutil.copy2(os.path.join(setup_files_liquid_dir, 'step_2.jou'), cfd2_dir)
        print(f"  Copied 'step_2.jou' to '{cfd2_dir}'")
    except Exception as e:
        print(f"  Error copying step_1 & step_2 jou files: {e}. Check setup_files/liquid/.")
        return

    # --- Step 7: Create 'remesh' folders ---
    print("\nStep 7: Creating 'remesh' folders...")
    cfd2_remesh_dir = os.path.join(cfd2_dir, 'remesh')
    try:
        if os.path.exists(cfd2_remesh_dir) and os.path.isdir(cfd2_remesh_dir):
            shutil.rmtree(cfd2_remesh_dir)
        os.makedirs(cfd2_remesh_dir, exist_ok=True)
        print(f"  Created '{cfd2_remesh_dir}'")
    except OSError as e:
        print(f"  Error handling remesh directory: {e}")
        return

    # --- Step 8: Handle Mesh Generation Strategy ---
    print("\nStep 8: Initialising mesh injection strategy...")
    if manual_mesh_path:
        if not os.path.exists(manual_mesh_path):
            print(f"  Error: The manual mesh file was not found at: {manual_mesh_path}")
            return
        try:
            dest_path = os.path.join(cfd2_remesh_dir, 'remesh.msh')
            shutil.copy2(manual_mesh_path, dest_path)
            print(f"  [MANUAL MODE] Copied '{manual_mesh_path}' to '{dest_path}'")
        except Exception as e:
            print(f"  Error copying manual mesh file: {e}")
            return
    else:
        remesh_files_2 = ['remesh.sh', 'identifyer.jou', 'generate_remesh.py']
        try:
            for f in remesh_files_2:
                shutil.copy2(os.path.join(setup_files_liquid_remesh_dir, f), cfd2_remesh_dir)
                print(f"  [AUTO MODE] Copied '{f}' to '{cfd2_remesh_dir}'")
        except Exception as e:
            print(f"  Error copying automated remesh files: {e}")
            return

    # --- Step 9: Prepare for Fluent step_1 run (copying latest timestep) ---
    print("\nStep 9: Preparing for Fluent step_1 run (copying latest timestep)...")

    def get_timestep_number(filepath):
        match = re.search(r'timestep(\d+)\.(cas|dat)\.h5$', os.path.basename(filepath))
        return int(match.group(1)) if match else -1

    # In manual mode, restore udf_thermal library and source from the archive if available
    if manual_mesh_path and max_run_num > 0:
        src_udf_dir = os.path.join(search_dir, 'udf_thermal')
        dst_udf_dir = os.path.join(cfd2_dir, 'udf_thermal')
        if os.path.exists(src_udf_dir):
            if os.path.exists(dst_udf_dir):
                shutil.rmtree(dst_udf_dir)
            shutil.copytree(src_udf_dir, dst_udf_dir)
            print(f"  [MANUAL MODE] Restored 'udf_thermal' library from {search_dir}")
        src_udf_c = os.path.join(search_dir, 'udf_thermal.c')
        if os.path.exists(src_udf_c):
            shutil.copy2(src_udf_c, os.path.join(cfd2_dir, 'udf_thermal.c'))
            print(f"  [MANUAL MODE] Restored 'udf_thermal.c' from {search_dir}")

    print(f"  Scanning for restart files in {search_dir}...")
    case_files = glob.glob(os.path.join(search_dir, 'case_timestep*.cas.h5'))
    data_files = glob.glob(os.path.join(search_dir, 'case_timestep*.dat.h5'))

    highest_case_file = sorted(case_files, key=get_timestep_number, reverse=True)[0] if case_files else None
    highest_data_file = sorted(data_files, key=get_timestep_number, reverse=True)[0] if data_files else None

    dest_case_restart = os.path.join(cfd2_dir, 'case_restart.cas.h5')
    dest_data_restart = os.path.join(cfd2_dir, 'case_restart.dat.h5')

    if highest_case_file and highest_data_file:
        try:
            shutil.copy2(highest_case_file, dest_case_restart)
            shutil.copy2(highest_data_file, dest_data_restart)
            print(f"  Copied latest timestep files ({os.path.basename(highest_case_file)}) to restart files in {cfd2_dir}")

            # Clean up older timestep files in the active directory if in automated mode
            if not manual_mesh_path:
                for f in sorted(case_files, key=get_timestep_number, reverse=True)[1:]:
                    os.remove(f)
                for f in sorted(data_files, key=get_timestep_number, reverse=True)[1:]:
                    os.remove(f)
        except Exception as e:
            print(f"  Error handling restart files: {e}")
    else:
        # Fallback: check for consolidated case_restart files in search_dir
        fallback_case = os.path.join(search_dir, 'case_restart.cas.h5')
        fallback_data = os.path.join(search_dir, 'case_restart.dat.h5')
        if os.path.exists(fallback_case) and os.path.exists(fallback_data):
            try:
                shutil.copy2(fallback_case, dest_case_restart)
                shutil.copy2(fallback_data, dest_data_restart)
                print(f"  Copied consolidated restart files from {search_dir} to {cfd2_dir}")
            except Exception as e:
                print(f"  Error copying fallback restart files: {e}")
        else:
            print(f"  Warning: Could not find previous case/data files in {search_dir}.")

    # --- Step 10: Run 'fluent' for 'step_1' ---
    print(f"\nStep 10: Running Fluent for 'step_1.jou'...")
    fluent_command_template = "ml -GAMBIT && ml ANSYS_CFD/2024R2 && fluent 2ddp -g -i {jou_file}"
    jou_file_name = "step_1.jou"
    full_jou_path = os.path.join(cfd2_dir, jou_file_name)
    command = fluent_command_template.format(jou_file=full_jou_path)

    try:
        result = subprocess.run(command, shell=True, cwd=cfd2_dir)
        if result.returncode != 0:
            print(f"  Fluent step_1 failed with code {result.returncode}.")
    except Exception as e:
        print(f"  Error running Fluent step_1: {e}")

    # --- Step 11: Post-step_1 cleanup ---
    print("\nStep 11: Cleaning up CFD_2 after step_1...")
    if is_first_run and not manual_mesh_path:
        stefan_folder = os.path.join(cfd2_dir, 'stefan')
        if os.path.exists(stefan_folder):
            shutil.rmtree(stefan_folder)
        if os.path.exists(os.path.join(cfd2_dir, 'stefan.c')):
            os.remove(os.path.join(cfd2_dir, 'stefan.c'))

    udf_thermal_folder = os.path.join(cfd2_dir, 'udf_thermal')
    if os.path.exists(udf_thermal_folder):
        shutil.rmtree(udf_thermal_folder)

    files_to_remove_post_step1 = ['v2024R2.jou', 'setup_fluent.log', 'setup_fluent_1.sh', 'udf_thermal.c']
    for f_name in files_to_remove_post_step1:
        f_path = os.path.join(cfd2_dir, f_name)
        if os.path.exists(f_path):
            os.remove(f_path)

    # --- Step 12: Execute Remeshing Script (If Automated) ---
    print("\nStep 12: Evaluating remesh script execution...")
    if manual_mesh_path:
        print("  Skipping remeshing script (Manual Mesh Used).")
    else:
        remesh_script_path = os.path.join(cfd2_remesh_dir, 'remesh.sh')
        if os.path.exists(remesh_script_path):
            print(f"  Executing automated remesh in {cfd2_remesh_dir}...")
            try:
                os.chmod(remesh_script_path, 0o755)
                result = subprocess.run(remesh_script_path, shell=True, cwd=cfd2_remesh_dir)
                if result.returncode != 0:
                    print(f"  remesh.sh failed with exit code {result.returncode}.")
            except Exception as e:
                print(f"  An error occurred while running remesh.sh: {e}")
        else:
            print(f"  Warning: remesh.sh not found in {cfd2_remesh_dir}. Skipping.")

    # --- Step 13: Run 'fluent' for 'step_2.jou' ---
    print("\nStep 13: Running Fluent for step_2.jou...")
    step2_jou_name = 'step_2.jou'
    full_jou_path = os.path.join(cfd2_dir, step2_jou_name)
    command = fluent_command_template.format(jou_file=full_jou_path)

    try:
        result = subprocess.run(command, shell=True, cwd=cfd2_dir)
        if result.returncode != 0:
            print(f"  Fluent step_2 failed with code {result.returncode}.")
    except Exception as e:
        print(f"  Error running Fluent step_2: {e}")

    # --- Step 14: Final cleanup ---
    print("\nStep 14: Performing final cleanup...")
    files_to_remove_final = ['*.trn', 'solver_load_cmd.log', 'case_timestep*', 'fluent.log', 'log']
    for pattern in files_to_remove_final:
        for file_path in glob.glob(os.path.join(cfd2_dir, pattern)):
            try:
                os.remove(file_path)
            except OSError:
                pass

    print("\n" + "=" * 50)
    mode = "MANUAL" if manual_mesh_path else "AUTOMATED"
    print(f"CFD Workflow Automation Completed ({mode} MODE)!")
    print("=" * 50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the CFD remeshing workflow.")
    parser.add_argument("--manual-mesh", type=str, default=None,
                        help="Path to a manually generated .msh file to override auto-remeshing.")

    args = parser.parse_args()
    run_cfd_workflow(args.manual_mesh)