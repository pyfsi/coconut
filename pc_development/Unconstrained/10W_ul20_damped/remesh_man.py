import os
import shutil
import glob
import re
import subprocess
import argparse
import sys


def run_cfd_workflow(manual_mesh_path):
    """
    Automates the CFD simulation workflow including file cleanup,
    directory management, and Fluent execution.
    """
    print(f"Starting CFD workflow automation script with manual mesh: {manual_mesh_path}")

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

    # --- Step 2: Determine run number ---
    print("\nStep 2: Determining run number and creating run folder...")
    existing_run_folders = glob.glob(os.path.join(current_dir, 'run_*'))
    max_run_num = 0
    for folder in existing_run_folders:
        match = re.match(r'run_(\d+)', os.path.basename(folder))
        if match:
            max_run_num = max(max_run_num, int(match.group(1)))

    previous_run_folder_name = f'run_{max_run_num}'
    previous_run_folder_path = os.path.join(current_dir, previous_run_folder_name)

    # --- Step 3: Remove local results files ---
    print("\nStep 3: Removing local results files...")
    files_to_remove_local = ['case_results.pickle', 'case_restart*']
    for pattern in files_to_remove_local:
        for file_path in glob.glob(os.path.join(current_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"  Removed: {file_path}")
            except OSError as e:
                print(f"  Error removing {file_path}: {e}")

    # --- Step 4 & 5: Copy 'step_1' and 'step_2' jou files ---
    print("\nStep 4 & 5: Copying Fluent journal files...")
    try:
        shutil.copy2(os.path.join(setup_files_liquid_dir, 'step_1.jou'), cfd2_dir)
        print(f"  Copied 'step_1.jou' to '{cfd2_dir}'")
        shutil.copy2(os.path.join(setup_files_liquid_dir, 'step_2.jou'), cfd2_dir)
        print(f"  Copied 'step_2.jou' to '{cfd2_dir}'")
    except Exception as e:
        print(f"  Error copying jou files: {e}")
        return

    # --- Step 6: Create 'remesh' folders ---
    print("\nStep 6: Creating 'remesh' folders...")
    cfd2_remesh_dir = os.path.join(cfd2_dir, 'remesh')

    try:
        if os.path.exists(cfd2_remesh_dir) and os.path.isdir(cfd2_remesh_dir):
            shutil.rmtree(cfd2_remesh_dir)
        os.makedirs(cfd2_remesh_dir, exist_ok=True)
        print(f"  Created '{cfd2_remesh_dir}'")
    except OSError as e:
        print(f"  Error handling remesh directory: {e}")
        return

    # --- Step 7: Copy 'remesh.msh' using the ARGUMENT ---
    print(f"\nStep 7: Copying manually remeshed mesh from argument: {manual_mesh_path}")

    if not os.path.exists(manual_mesh_path):
        print(f"  Error: The manual mesh file was not found at: {manual_mesh_path}")
        return

    try:
        # Copy the file provided in arguments to the destination as 'remesh.msh'
        dest_path = os.path.join(cfd2_remesh_dir, 'remesh.msh')
        shutil.copy2(manual_mesh_path, dest_path)
        print(f"  Copied '{manual_mesh_path}' to '{dest_path}'")
    except Exception as e:
        print(f"  Error copying manual mesh file: {e}")
        return

    # --- Step 8: Prepare for Fluent step_1 run ---
    print("\nStep 8: Preparing for Fluent step_1 run (copying latest timestep)...")

    def get_timestep_number(filepath):
        match = re.search(r'timestep(\d+)\.(cas|dat)\.h5$', os.path.basename(filepath))
        return int(match.group(1)) if match else -1

    prev_cfd_dir = os.path.join(previous_run_folder_path, 'CFD_2')
    case_files = glob.glob(os.path.join(prev_cfd_dir, 'case_timestep*.cas.h5'))
    data_files = glob.glob(os.path.join(prev_cfd_dir, 'case_timestep*.dat.h5'))

    highest_case_file = sorted(case_files, key=get_timestep_number, reverse=True)[0] if case_files else None
    highest_data_file = sorted(data_files, key=get_timestep_number, reverse=True)[0] if data_files else None

    if highest_case_file and highest_data_file:
        # Fixed: cfd_dir_2 was undefined, changed to cfd2_dir
        dest_case_restart = os.path.join(cfd2_dir, 'case_restart.cas.h5')
        dest_data_restart = os.path.join(cfd2_dir, 'case_restart.dat.h5')
        try:
            shutil.copy2(highest_case_file, dest_case_restart)
            shutil.copy2(highest_data_file, dest_data_restart)
            print(f"  Copied latest timestep files to restart files in {cfd2_dir}")
        except Exception as e:
            print(f"  Error copying restart files: {e}")
    else:
        print(f"  Warning: Could not find previous case/data files in {prev_cfd_dir}.")

    # --- Step 9: Run 'fluent' for 'step_1' ---
    print(f"\nStep 9: Running Fluent for 'step_1.jou'...")
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

    # --- Step 10: Post-step_1 cleanup ---
    print("\nStep 10: Cleaning up CFD_2 after step_1...")
    udf_thermal_folder = os.path.join(cfd2_dir, 'udf_thermal')
    if os.path.exists(udf_thermal_folder):
        shutil.rmtree(udf_thermal_folder)

    files_to_remove_post_step1 = ['v2024R2.jou', 'setup_fluent.log', 'setup_fluent_1.sh', 'udf_thermal.c']
    for f_name in files_to_remove_post_step1:
        f_path = os.path.join(cfd2_dir, f_name)
        if os.path.exists(f_path):
            os.remove(f_path)

    # --- Step 11: Skip remeshing script ---
    print("\nStep 11: Skipping remeshing script (Manual Mesh Used).")

    # --- Step 12: Run 'fluent' for 'step_2.jou' ---
    print("\nStep 12: Running Fluent for step_2.jou...")
    step2_jou_name = 'step_2.jou'
    full_jou_path = os.path.join(cfd2_dir, step2_jou_name)
    command = fluent_command_template.format(jou_file=full_jou_path)

    try:
        result = subprocess.run(command, shell=True, cwd=cfd2_dir)
        if result.returncode != 0:
            print(f"  Fluent step_2 failed with code {result.returncode}.")
    except Exception as e:
        print(f"  Error running Fluent step_2: {e}")

    # --- Step 13: Final cleanup ---
    print("\nStep 13: Performing final cleanup...")
    files_to_remove_final = ['*.trn', 'solver_load_cmd.log', 'case_timestep*', 'fluent.log', 'log']
    for pattern in files_to_remove_final:
        for file_path in glob.glob(os.path.join(cfd2_dir, pattern)):
            try:
                os.remove(file_path)
            except OSError:
                pass

    print("\n" + "=" * 50)
    print("CFD Workflow Automation Completed (Manual Remesh)!")
    print("=" * 50)


if __name__ == "__main__":
    # Parsing arguments
    parser = argparse.ArgumentParser(description="Run CFD workflow with manual mesh.")
    parser.add_argument("manual_mesh", help="Path to the manual remesh file (.msh)")
    args = parser.parse_args()

    run_cfd_workflow(args.manual_mesh)