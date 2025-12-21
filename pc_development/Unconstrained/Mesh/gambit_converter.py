import csv
import os

# --- Configuration ---
# Your original CSV file
input_file_name = 'sorted_interface.csv'
# ---------------------

# Generate the new output filename
base_name, extension = os.path.splitext(input_file_name)
output_file_name = f"{base_name}_gambit{extension}"

print(f"Starting conversion...")
print(f"  Input:  {input_file_name}")
print(f"  Output: {output_file_name}")

try:
    with open(input_file_name, mode='r', newline='') as infile, \
            open(output_file_name, mode='w', newline='') as outfile:

        csv_reader = csv.reader(infile)

        # 1. Skip the header row
        try:
            next(csv_reader)
        except StopIteration:
            print(f"Error: Input file '{input_file_name}' appears to be empty.")
            exit()  # Stop if there's no data

        # 2. Read data rows, reformat, and write to the new file
        count = 0
        for row in csv_reader:
            if len(row) >= 2:  # Check for at least two columns
                x = row[0]  # First column
                y = row[1]  # Second column
                z = "0.0"  # Hardcoded z-value for 2D geometry

                # Write in the space-delimited "x y z" format
                outfile.write(f"{x} {y} {z}\n")
                count += 1
            else:
                print(f"Skipping malformed row: {row}")

    print(f"\nSuccessfully converted {count} vertices.")

except FileNotFoundError:
    print(f"\nError: The file '{input_file_name}' was not found.")
    print("Please make sure the script is in the same directory as your CSV file.")
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}")