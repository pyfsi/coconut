import os
import shutil
import glob
import re
import subprocess

def run_cfd_workflow():
    """
    Automates the CFD simulation workflow including file cleanup,
    directory management, and Fluent execution.
    """
    print("Starting CFD workflow automation script...")

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
        '*.coco',
        '*.trn',
        'cleanup-*.sh'
        'displacement_timestep*',
        'faces_timestep*',
        'heat_flux_timestep*',
        'nodes*'
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

    # --- Step 2: Determine run number and create 'run_x' folder ---
    print("\nStep 2: Determining run number and creating run folder...")
    existing_run_folders = glob.glob(os.path.join(current_dir, 'run_*'))
    max_run_num = 0
    for folder in existing_run_folders:
        match = re.match(r'run_(\d+)', os.path.basename(folder))
        if match:
            max_run_num = max(max_run_num, int(match.group(1)))

    next_run_num = max_run_num + 1
    new_run_folder_name = f'run_{next_run_num}'
    new_run_folder_path = os.path.join(current_dir, new_run_folder_name)

    is_first_run = (next_run_num == 1)

    try:
        os.makedirs(new_run_folder_path, exist_ok=True)
        print(f"  Created run folder: {new_run_folder_path}")
    except OSError as e:
        print(f"  Error creating folder {new_run_folder_path}: {e}")
        return

    # --- Step 3: Copy CFD_1, CFD_2, and results files to the new 'run_x' folder ---
    print(f"\nStep 3: Copying CFD folders and results to {new_run_folder_name}...")
    files_to_copy_to_run = [
        'case_results.pickle',
        'case_restart*'
    ]

    # Copy CFD_1 and CFD_2 folders
    for src_cfd_dir in [cfd1_dir, cfd2_dir]:
        dest_cfd_dir = os.path.join(new_run_folder_path, os.path.basename(src_cfd_dir))
        try:
            # shutil.copytree requires destination to not exist, so we delete if it does
            if os.path.exists(dest_cfd_dir):
                shutil.rmtree(dest_cfd_dir)
            shutil.copytree(src_cfd_dir, dest_cfd_dir)
            print(f"  Copied '{os.path.basename(src_cfd_dir)}' to '{new_run_folder_name}'")
        except shutil.Error as e:
            print(f"  Error copying directory {src_cfd_dir} to {dest_cfd_dir}: {e}")

    # Copy specified files from current_dir to new_run_folder_path
    for pattern in files_to_copy_to_run:
        for file_path in glob.glob(os.path.join(current_dir, pattern)):
            try:
                shutil.copy2(file_path, new_run_folder_path)
                print(f"  Copied '{os.path.basename(file_path)}' to '{new_run_folder_name}'")
            except FileNotFoundError:
                print(f"  Warning: File matching '{pattern}' not found in current directory. Skipping.")
            except Exception as e:
                print(f"  Error copying {file_path}: {e}")

    # --- Step 4: Remove 'case_results.pickle' and 'case_restart*' from current directory ---
    print("\nStep 4: Removing local results files...")
    files_to_remove_local = [
        'case_results.pickle',
        'case_restart*'
    ]
    for pattern in files_to_remove_local:
        for file_path in glob.glob(os.path.join(current_dir, pattern)):
            try:
                os.remove(file_path)
                print(f"  Removed: {file_path}")
            except OSError as e:
                print(f"  Error removing {file_path}: {e}")

    # --- Step 5 & 6: Copy 'step_1' and 'step_2' jou files ---
    print("\nStep 5 & 6: Copying Fluent journal files...")
    step1_jou_name = 'step_1_firstrun.jou' if is_first_run else 'step_1_nextrun.jou'

    # Copy step_1_*.jou
    try:
        shutil.copy2(os.path.join(setup_files_solid_dir, step1_jou_name), cfd1_dir)
        print(f"  Copied '{step1_jou_name}' to '{cfd1_dir}'")
        shutil.copy2(os.path.join(setup_files_liquid_dir, step1_jou_name), cfd2_dir)
        print(f"  Copied '{step1_jou_name}' to '{cfd2_dir}'")
    except FileNotFoundError as e:
        print(f"  Error: {e}. Make sure '{step1_jou_name}' exists in 'setup_files/solid/' and 'setup_files/liquid/'.")
        return
    except Exception as e:
        print(f"  Error copying step_1 jou files: {e}")
        return

    # Copy step_2.jou
    try:
        shutil.copy2(os.path.join(setup_files_solid_dir, 'step_2.jou'), cfd1_dir)
        print(f"  Copied 'step_2.jou' to '{cfd1_dir}'")
        shutil.copy2(os.path.join(setup_files_liquid_dir, 'step_2.jou'), cfd2_dir)
        print(f"  Copied 'step_2.jou' to '{cfd2_dir}'")
    except FileNotFoundError as e:
        print(f"  Error: {e}. Make sure 'step_2.jou' exists in 'setup_files/solid/' and 'setup_files/liquid/'.")
        return
    except Exception as e:
        print(f"  Error copying step_2 jou files: {e}")
        return

    # --- Step 7: Create 'remesh' folders ---
    print("\nStep 7: Creating 'remesh' folders...")
    cfd1_remesh_dir = os.path.join(cfd1_dir, 'remesh')
    cfd2_remesh_dir = os.path.join(cfd2_dir, 'remesh')

    for remesh_dir in [cfd1_remesh_dir, cfd2_remesh_dir]:
        try:
            # Check if remesh directory exists and delete it
            if os.path.exists(remesh_dir) and os.path.isdir(remesh_dir):
                print(f"  Existing '{remesh_dir}' found. Deleting...")
                shutil.rmtree(remesh_dir)
                print(f"  Deleted existing '{remesh_dir}'")

            os.makedirs(remesh_dir, exist_ok=True)
            print(f"  Created '{remesh_dir}'")
        except OSError as e:
            print(f"  Error creating remesh directory {remesh_dir}: {e}")
            return
        except shutil.Error as e:
            print(f"  Error deleting existing remesh directory {remesh_dir}: {e}")
            return

    # --- Step 8: Copy 'remesh.sh' and 'remesh.jou' ---
    print("\nStep 8: Copying remesh files...")
    remesh_files = ['remesh.sh', 'remesh.jou']

    try:
        for f in remesh_files:
            shutil.copy2(os.path.join(setup_files_solid_remesh_dir, f), cfd1_remesh_dir)
            print(f"  Copied '{f}' to '{cfd1_remesh_dir}'")
            shutil.copy2(os.path.join(setup_files_liquid_remesh_dir, f), cfd2_remesh_dir)
            print(f"  Copied '{f}' to '{cfd2_remesh_dir}'")
    except FileNotFoundError as e:
        print(f"  Error: {e}. Make sure 'remesh.sh' and 'remesh.jou' exist in 'setup_files/solid/create_mesh/' and 'setup_files/liquid/create_mesh/'.")
        return
    except Exception as e:
        print(f"  Error copying remesh files: {e}")
        return

    # --- Step 9: Prepare for Fluent step_1 run (copying latest timestep to restart files) ---
    print("\nStep 9: Preparing for Fluent step_1 run (copying latest timestep to restart files)...")

    # Function to extract timestep number
    def get_timestep_number(filepath):
        match = re.search(r'timestep(\d+)\.(cas|dat)\.h5$', os.path.basename(filepath))
        return int(match.group(1)) if match else -1

    for cfd_dir in [cfd1_dir, cfd2_dir]:
        print(f"  Processing timestep files in {cfd_dir}...")
        case_files = glob.glob(os.path.join(cfd_dir, 'case_timestep*.cas.h5'))
        data_files = glob.glob(os.path.join(cfd_dir, 'case_timestep*.dat.h5'))

        highest_case_file = None
        highest_data_file = None

        # Process .cas.h5 files
        if case_files:
            case_files_sorted = sorted(case_files, key=get_timestep_number, reverse=True)
            highest_case_file = case_files_sorted[0]
            print(f"    Found latest case file: {highest_case_file}")
            # Remove older case files
            for f in case_files_sorted[1:]:
                try:
                    os.remove(f)
                    print(f"      Removed old case file: {f}")
                except OSError as e:
                    print(f"      Error removing {f}: {e}")
        else:
            print(f"    No 'case_timestep*.cas.h5' files found in {cfd_dir}.")

        # Process .dat.h5 files
        if data_files:
            data_files_sorted = sorted(data_files, key=get_timestep_number, reverse=True)
            highest_data_file = data_files_sorted[0]
            print(f"    Found latest data file: {highest_data_file}")
            # Remove older data files
            for f in data_files_sorted[1:]:
                try:
                    os.remove(f)
                    print(f"      Removed old data file: {f}")
                except OSError as e:
                    print(f"      Error removing {f}: {e}")
        else:
            print(f"    No 'case_timestep*.dat.h5' files found in {cfd_dir}.")

        # Copy highest timestep files to case_restart.cas.h5 and case_restart.dat.h5
        if highest_case_file and highest_data_file:
            dest_case_restart = os.path.join(cfd_dir, 'case_restart.cas.h5')
            dest_data_restart = os.path.join(cfd_dir, 'case_restart.dat.h5')
            try:
                shutil.copy2(highest_case_file, dest_case_restart)
                print(
                    f"    Copied '{os.path.basename(highest_case_file)}' to '{os.path.basename(dest_case_restart)}'")
                shutil.copy2(highest_data_file, dest_data_restart)
                print(
                    f"    Copied '{os.path.basename(highest_data_file)}' to '{os.path.basename(dest_data_restart)}'")
            except Exception as e:
                print(f"    Error copying latest timestep files to restart files in {cfd_dir}: {e}")
        else:  # This block will now always execute if files are missing, even for first run
            print(
                f"    Warning: Could not find both latest case and data files to create restart files in {cfd_dir}. Ensure previous run generated them.")

    # --- Step 10: Run 'fluent' for 'step_1' ---
    print(f"\nStep 10: Running Fluent for {step1_jou_name}...")
    fluent_command_template = "ml -GAMBIT && ml ANSYS_CFD/2024R2 && fluent 2ddp -g -i {jou_file}"

    for cfd_dir, jou_file_name in [(cfd1_dir, step1_jou_name), (cfd2_dir, step1_jou_name)]:
        full_jou_path = os.path.join(cfd_dir, jou_file_name)
        command = fluent_command_template.format(jou_file=full_jou_path)
        print(f"  Executing in {cfd_dir}: {command}")
        try:
            # Use subprocess.run for better control and error handling
            # cwd sets the current working directory for the command
            result = subprocess.run(command, shell=True, cwd=cfd_dir)
            if result.returncode == 0:
                print(f"  Fluent run for {jou_file_name} in {cfd_dir} completed successfully.")
                # print("STDOUT:\n", result.stdout) # Uncomment to see Fluent output
            else:
                print(f"  Fluent run for {jou_file_name} in {cfd_dir} failed with exit code {result.returncode}.")
                print("  STDERR:\n", result.stderr)
        except FileNotFoundError:
            print(f"  Error: 'fluent' command not found. Ensure Fluent is installed and in your PATH.")
        except Exception as e:
            print(f"  An error occurred while running Fluent for {jou_file_name} in {cfd_dir}: {e}")

    # --- Step 11: Post-step_1 cleanup in CFD_1 and CFD_2 (excluding timestep files, now handled in 8.5) ---
    print("\nStep 11: Cleaning up CFD_1 and CFD_2 after step_1 Fluent run (excluding timestep files)...")
    for cfd_dir in [cfd1_dir, cfd2_dir]:
        print(f"  Cleaning in {cfd_dir}...")
        # Remove 'stefan' folder and file if first run
        if is_first_run:
            stefan_folder = os.path.join(cfd_dir, 'stefan')
            stefan_file = os.path.join(cfd_dir, 'stefan.c')
            if os.path.exists(stefan_folder):
                try:
                    shutil.rmtree(stefan_folder)
                    print(f"    Removed folder: {stefan_folder}")
                except OSError as e:
                    print(f"    Error removing {stefan_folder}: {e}")
            if os.path.exists(stefan_file):
                try:
                    os.remove(stefan_file)
                    print(f"    Removed file: {stefan_file}")
                except OSError as e:
                    print(f"    Error removing {stefan_file}: {e}")

        # Remove 'udf_thermal.c' and folder 'udf_thermal'
        udf_thermal_folder = os.path.join(cfd_dir, 'udf_thermal')
        udf_thermal_file = os.path.join(cfd_dir, 'udf_thermal.c')
        if os.path.exists(udf_thermal_folder):
            try:
                shutil.rmtree(udf_thermal_folder)
                print(f"    Removed folder: {udf_thermal_folder}")
            except OSError as e:
                print(f"    Error removing {udf_thermal_folder}: {e}")
        if os.path.exists(udf_thermal_file):
            try:
                os.remove(udf_thermal_file)
                print(f"    Removed file: {udf_thermal_file}")
            except OSError as e:
                print(f"    Error removing {udf_thermal_file}: {e}")

        # Remove 'v2024R2.jou', 'setup_fluent.log', 'setup_fluent_1.sh'
        files_to_remove_post_step1 = [
            'v2024R2.jou',
            'setup_fluent.log',
            'setup_fluent_1.sh'
        ]
        for f_name in files_to_remove_post_step1:
            file_path = os.path.join(cfd_dir, f_name)
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                    print(f"    Removed: {file_path}")
                except OSError as e:
                    print(f"    Error removing {file_path}: {e}")

        # The timestep cleanup logic was moved to Step 8.5

    # --- Step 12: Run remesh.sh shell script ---
    print("\nStep 12: Running remesh.sh scripts...")
    for cfd_remesh_dir in [cfd1_remesh_dir, cfd2_remesh_dir]:
        remesh_script_path = os.path.join(cfd_remesh_dir, 'remesh.sh')
        if os.path.exists(remesh_script_path):
            print(f"  Executing in {cfd_remesh_dir}: {remesh_script_path}")
            try:
                # Ensure the script is executable
                os.chmod(remesh_script_path, 0o755)
                result = subprocess.run(remesh_script_path, shell=True, cwd=cfd_remesh_dir)
                if result.returncode == 0:
                    print(f"  remesh.sh in {cfd_remesh_dir} completed successfully.")
                    # print("STDOUT:\n", result.stdout) # Uncomment to see remesh output
                else:
                    print(f"  remesh.sh in {cfd_remesh_dir} failed with exit code {result.returncode}.")
                    print("  STDERR:\n", result.stderr)
            except Exception as e:
                print(f"  An error occurred while running remesh.sh in {cfd_remesh_dir}: {e}")
        else:
            print(f"  Warning: remesh.sh not found in {cfd_remesh_dir}. Skipping.")

    # --- Step 13: Run 'fluent 2ddp -g -i step_2.jou' ---
    print("\nStep 13: Running Fluent for step_2.jou...")
    step2_jou_name = 'step_2.jou'
    for cfd_dir in [cfd1_dir, cfd2_dir]:
        full_jou_path = os.path.join(cfd_dir, step2_jou_name)
        command = fluent_command_template.format(jou_file=full_jou_path)
        print(f"  Executing in {cfd_dir}: {command}")
        try:
            result = subprocess.run(command, shell=True, cwd=cfd_dir)
            if result.returncode == 0:
                print(f"  Fluent run for {step2_jou_name} in {cfd_dir} completed successfully.")
                # print("STDOUT:\n", result.stdout) # Uncomment to see Fluent output
            else:
                print(f"  Fluent run for {step2_jou_name} in {cfd_dir} failed with exit code {result.returncode}.")
                print("  STDERR:\n", result.stderr)
        except FileNotFoundError:
            print(f"  Error: 'fluent' command not found. Ensure Fluent is installed and in your PATH.")
        except Exception as e:
            print(f"  An error occurred while running Fluent for {step2_jou_name} in {cfd_dir}: {e}")

    # --- Step 14: Final cleanup in CFD_1 and CFD_2 ---
    print("\nStep 14: Performing final cleanup in CFD_1 and CFD_2...")
    files_to_remove_final = [
        '*.trn',
        'solver_load_cmd.log'
    ]
    for cfd_dir in [cfd1_dir, cfd2_dir]:
        print(f"  Cleaning in {cfd_dir}...")
        for pattern in files_to_remove_final:
            for file_path in glob.glob(os.path.join(cfd_dir, pattern)):
                try:
                    os.remove(file_path)
                    print(f"    Removed: {file_path}")
                except OSError as e:
                    print(f"    Error removing {file_path}: {e}")

    # --- Step 15: Print finishing message ---
    print("\n" + "="*50)
    print("CFD Workflow Automation Completed!")
    print(f"New run created: {new_run_folder_name}")
    print("\nREMINDER:")
    print("- Update 'report-file.out' with the latest results.")
    print("- Update the naming of the Fluent case and data files if necessary.")
    print("- Update the 'parameters.json' file for the next simulation.")
    print("="*50)

if __name__ == "__main__":
    run_cfd_workflow()
