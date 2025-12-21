import csv
import os

# --- Configuration ---
# Your original CSV file
input_file_name = 'sorted_interface.csv'
# The constant group ID you want for this curve
group_id = 1
# ---------------------

# Generate the new output filename for DesignModeler
base_name, extension = os.path.splitext(input_file_name)
# Use .txt extension, as this is a specific text-based format, not a CSV
output_file_name = f"{base_name}_designmodeler.txt"

print(f"Starting conversion for DesignModeler (3D Curve format)...")
print(f"  Input:  {input_file_name}")
print(f"  Output: {output_file_name}")

try:
    with open(input_file_name, mode='r', newline='') as infile, \
            open(output_file_name, mode='w', newline='') as outfile:

        csv_reader = csv.reader(infile)

        # 1. Skip the header row from the input CSV
        try:
            next(csv_reader)
        except StopIteration:
            print(f"Error: Input file '{input_file_name}' appears to be empty.")
            exit()  # Stop if there's no data

        # 2. Write data rows in the new format
        point_index = 1
        for row in csv_reader:
            if len(row) >= 2:
                x = row[0].strip()  # Get x from 1st column
                y = row[1].strip()  # Get y from 2nd column
                z = "0"  # Z-coordinate is 0, as per your example

                # Write in the format: [GroupID] [PointID] [X] [Y] [Z]
                outfile.write(f"{group_id} {point_index} {x} {y} {z}\n")
                point_index += 1
            else:
                print(f"Skipping malformed row: {row}")

        # 3. Write the terminator line for the group
        outfile.write(f"{group_id} 0\n")

        # We subtract 1 from point_index because it was incremented one last time
        total_points = point_index - 1
        print(f"\nSuccessfully converted {total_points} vertices into Group {group_id}.")

except FileNotFoundError:
    print(f"\nError: The file '{input_file_name}' was not found.")
    print("Please make sure the script is in the same directory as your CSV file.")
except Exception as e:
    print(f"\nAn unexpected error occurred: {e}")