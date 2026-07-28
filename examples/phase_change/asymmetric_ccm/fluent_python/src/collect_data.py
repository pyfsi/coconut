import os
import shutil
import pandas as pd
import numpy as np

# === CONFIGURATION ===
CFD_DIR = "CFD_2"
POST_PROC_DIR = "post_processing"  # Updated to match your log
REPORT_FILE = "report-file.out"
RB_FILE = "rigid-body-report-file.out"


def ensure_post_proc_dir():
    if not os.path.exists(POST_PROC_DIR):
        try:
            os.makedirs(POST_PROC_DIR)
            print(f"[CollectData] Created directory: {POST_PROC_DIR}")
        except OSError as e:
            print(f"[CollectData] Error creating directory {POST_PROC_DIR}: {e}")


def process_report_file():
    """
    Handles report-file.out:
    - SKIPS first 3 header lines to correctly detect columns.
    - DISCARDS Row 0 (Overlap).
    - Appends remainder with accumulated Time Step and Flow Time.
    """
    local_path = os.path.join(CFD_DIR, REPORT_FILE)
    global_path = os.path.join(POST_PROC_DIR, REPORT_FILE)

    if not os.path.exists(local_path):
        print(f"[CollectData] Warning: Local file {local_path} not found.")
        return

    try:
        # 1. Read Local File
        # FIX: skiprows=3 tells pandas to ignore the 3 header lines and infer columns from data
        df_local = pd.read_csv(local_path, sep=r'\s+', header=None, skiprows=3, on_bad_lines='skip')

        # Ensure numeric
        df_local = df_local[pd.to_numeric(df_local[0], errors='coerce').notnull()].astype(float)

        if df_local.empty:
            print(f"[CollectData] CRITICAL: Local {REPORT_FILE} is empty after filtering.")
            return

        if not os.path.exists(global_path):
            # --- FIRST RUN: INITIALIZE ---
            # Copy the full original file (headers + data) to preserve format
            shutil.copy2(local_path, global_path)
            print(f"[CollectData] Initialized {REPORT_FILE}")
        else:
            # --- APPEND MODE ---
            # 1. Get last values from Global
            # Read global with same skiprows logic to interpret columns correctly
            df_global = pd.read_csv(global_path, sep=r'\s+', header=None, skiprows=3, on_bad_lines='skip')
            df_global = df_global[pd.to_numeric(df_global[0], errors='coerce').notnull()].astype(float)

            if df_global.empty:
                # Fallback if global exists but is broken/empty
                shutil.copy2(local_path, global_path)
                return

            last_row = df_global.iloc[-1]
            last_step = last_row.iloc[0]  # Col 0: Time Step
            last_time = last_row.iloc[-1]  # Last Col: Flow Time

            # 2. Process Local: DISCARD Row 0 (The overlap)
            df_local = df_local[df_local[0] > 0]

            if not df_local.empty:
                # 3. Accumulate (Simple Offset)
                df_local[0] = df_local[0] + last_step
                df_local.iloc[:, -1] = df_local.iloc[:, -1] + last_time

                # 4. Append
                with open(global_path, 'a') as f:
                    df_local.to_csv(f, sep=' ', index=False, header=False)
                print(f"[CollectData] Appended {len(df_local)} rows to {REPORT_FILE}")
            else:
                print(f"[CollectData] No new data in {REPORT_FILE} (after removing step 0)")

    except Exception as e:
        print(f"[CollectData] Error processing {REPORT_FILE}: {e}")
        # import traceback
        # traceback.print_exc()


def format_rb_line(row):
    """Custom formatter for rigid body lines."""
    line_str = ""
    for i, val in enumerate(row):
        val_str = "{:.5e}".format(val)
        if i == 0:
            line_str += " " + val_str
        else:
            if val >= 0:
                line_str += "   " + val_str
            else:
                line_str += "  " + val_str
    return line_str


def process_rigid_body_file():
    """
    Handles rigid-body-report-file.out:
    - KEEPS Row 0.
    - Universal Time Offset.
    """
    local_path = os.path.join(CFD_DIR, RB_FILE)
    global_path = os.path.join(POST_PROC_DIR, RB_FILE)

    if not os.path.exists(local_path):
        print(f"[CollectData] Warning: Local file {local_path} not found.")
        return

    try:
        # 1. Read Local File
        df_local = pd.read_csv(local_path, sep=r'\s+', comment='#', header=None).astype(float)

        if not os.path.exists(global_path):
            # --- FIRST RUN ---
            shutil.copy2(local_path, global_path)
            print(f"[CollectData] Initialized {RB_FILE}")
        else:
            # --- APPEND MODE ---
            df_global = pd.read_csv(global_path, sep=r'\s+', comment='#', header=None).astype(float)

            if df_global.empty:
                shutil.copy2(local_path, global_path)
                return

            last_row = df_global.iloc[-1]
            last_time = last_row.iloc[0]  # Col 0: Time
            last_theta = last_row.iloc[5]  # Col 5: Theta

            # 2. Calculate Time Step (dt)
            dt = 0.0
            if len(df_global) >= 2:
                dt = df_global.iloc[-1, 0] - df_global.iloc[-2, 0]
            elif len(df_local) >= 2:
                dt = df_local.iloc[1, 0] - df_local.iloc[0, 0]

            if dt <= 0: dt = 0.0001

            # 3. Process Local: KEEP Row 0
            if not df_local.empty:
                # 4. Universal Offset Logic
                target_start_time = last_time + dt
                current_start_time = df_local.iloc[0, 0]
                time_offset = target_start_time - current_start_time

                # Apply Offsets
                df_local[0] = df_local[0] + time_offset
                df_local[5] = df_local[5] + last_theta

                # 5. Append with Custom Formatting
                with open(global_path, 'a') as f:
                    for index, row in df_local.iterrows():
                        formatted_line = format_rb_line(row)
                        f.write(formatted_line + "\n")

                print(f"[CollectData] Appended {len(df_local)} rows to {RB_FILE} (dt={dt:.2e})")

    except Exception as e:
        print(f"[CollectData] Error processing {RB_FILE}: {e}")


if __name__ == "__main__":
    print("--- Starting Data Collection ---")
    ensure_post_proc_dir()
    process_report_file()
    process_rigid_body_file()
    print("--- Data Collection Complete ---")