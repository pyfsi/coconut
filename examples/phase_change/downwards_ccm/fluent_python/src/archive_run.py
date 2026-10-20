import os
import shutil
import glob
import json
import re


def load_state():
    if os.path.exists('state.json'):
        with open('state.json', 'r') as file:
            return json.load(file)
    return {"current_run": 0}


def archive_final_run():
    current_dir = os.getcwd()
    state = load_state()
    run_num = state.get('current_run', 0)

    # Fallback to scanning existing folders if state is uninitialised
    if run_num == 0:
        existing_run_folders = glob.glob(os.path.join(current_dir, 'run_*'))
        max_run_num = 0
        for folder in existing_run_folders:
            match = re.match(r'run_(\d+)', os.path.basename(folder))
            if match:
                max_run_num = max(max_run_num, int(match.group(1)))
        run_num = max_run_num + 1

    run_folder_name = f'run_{run_num}'
    run_folder_path = os.path.join(current_dir, run_folder_name)

    print(f"\n[Archive] Archiving final run state to '{run_folder_name}'...")

    if os.path.exists(run_folder_path):
        print(f"[Archive] Warning: '{run_folder_name}' already exists. Updating contents...")
    else:
        os.makedirs(run_folder_path, exist_ok=True)
        print(f"[Archive] Created archive directory: {run_folder_path}")

    # 1. Copy CFD_1 and CFD_2
    for cfd_name in ['CFD_1', 'CFD_2']:
        src_path = os.path.join(current_dir, cfd_name)
        dest_path = os.path.join(run_folder_path, cfd_name)
        if os.path.exists(src_path):
            if os.path.exists(dest_path):
                shutil.rmtree(dest_path)
            shutil.copytree(src_path, dest_path)
            print(f"[Archive] Copied '{cfd_name}' to '{run_folder_name}'")
        else:
            print(f"[Archive] Warning: '{cfd_name}' directory not found at {src_path}")

    # 2. Copy restart and results pickle files
    files_to_copy = ['case_results.pickle', 'case_restart*', 'restart_rb_*.pickle']
    copied_count = 0
    for pattern in files_to_copy:
        for file_path in glob.glob(os.path.join(current_dir, pattern)):
            try:
                shutil.copy2(file_path, run_folder_path)
                print(f"[Archive] Copied '{os.path.basename(file_path)}' to '{run_folder_name}'")
                copied_count += 1
            except Exception as e:
                print(f"[Archive] Error copying {file_path}: {e}")

    print(f"[Archive] Final run archive completed: {run_folder_name} (Files copied: {copied_count})\n")


if __name__ == "__main__":
    archive_final_run()