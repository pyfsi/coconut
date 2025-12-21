import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
# Threshold factor for point simplification. If a segment is shorter than
# this factor times the average segment length, a node is removed.
# SIMPLIFICATION_THRESHOLD_FACTOR = 0.5
SIMPLIFICATION_THRESHOLD_FACTOR = 0

# Set a maximum Y-value for clamping points.
# All points with a y-coordinate > MAX_Y_VALUE will have their y-value set to MAX_Y_VALUE.
# Set to 'None' to disable clamping.
# MAX_Y_VALUE = 0.05076
MAX_Y_VALUE = None

# --- Helper function for distance calculation ---
def calculate_distance(p1, p2):
    """Calculates the Euclidean distance between two points."""
    return np.sqrt(np.sum((p1 - p2) ** 2))


# --- Main script ---
# This script now reads from an 'interface.xy' file.
file_path = 'interface_100.xy'

try:
    # --- Read and Parse the .xy file ---
    points_list = []
    with open(file_path, 'r') as f:
        for line in f:
            # Clean up the line
            line = line.strip()
            # Skip header lines (which start with '(') or empty lines
            if not line or line.startswith('('):
                continue

            try:
                # Split by whitespace and convert to float
                parts = line.split()
                if len(parts) == 2:
                    x = float(parts[0])
                    y = float(parts[1])
                    points_list.append([x, y])
            except (ValueError, IndexError):
                # This will skip any lines that aren't two valid numbers
                print(f"Skipping malformed line: {line}")
                continue

    # --- Clamp points based on MAX_Y_VALUE ---
    if MAX_Y_VALUE is not None:
        clamped_points_list = []
        clamped_count = 0
        for p in points_list:
            if p[1] > MAX_Y_VALUE:
                clamped_points_list.append([p[0], MAX_Y_VALUE])
                clamped_count += 1
            else:
                clamped_points_list.append(p)
        if clamped_count > 0:
            print(f"Clamped {clamped_count} points to y = {MAX_Y_VALUE}")
        points_to_convert = clamped_points_list
    else:
        points_to_convert = points_list

    # Convert the list of points into a NumPy array for efficient calculations.
    points = np.array(points_to_convert)

    if points.shape[0] == 0:
        raise ValueError("No valid coordinate data found in the input file.")

    # --- Nearest Neighbor Sorting ---

    n_points = len(points)
    if n_points > 1:
        # A list to store the path of indices.
        path = [0]
        # A boolean array to keep track of visited points.
        unvisited = np.ones(n_points, dtype=bool)
        unvisited[0] = False

        for i in range(n_points - 1):
            last_point_index = path[-1]
            last_point = points[last_point_index]

            # Calculate squared Euclidean distance to all unvisited points.
            # Using squared distance is faster as it avoids the square root.
            deltas = points[unvisited] - last_point
            dist_sq = np.sum(deltas ** 2, axis=1)

            # Find the index of the nearest point among the unvisited ones.
            nearest_local_index = np.argmin(dist_sq)

            # Get the original index of that nearest point.
            unvisited_indices = np.where(unvisited)[0]
            nearest_global_index = unvisited_indices[nearest_local_index]

            # Add the nearest point to the path and mark it as visited.
            path.append(nearest_global_index)
            unvisited[nearest_global_index] = False

        # Reorder the original points based on the nearest neighbor path.
        sorted_points = points[path]
    else:
        sorted_points = points

    # --- Point Simplification ---
    if len(sorted_points) > 3:  # Need at least a quadrilateral to simplify
        # Convert to a list for easier element removal
        simplified_points_list = list(sorted_points)

        while True:
            point_removed_this_pass = False
            num_points = len(simplified_points_list)

            if num_points <= 3:
                break

            # Calculate segment lengths for the current set of points
            lengths = [calculate_distance(simplified_points_list[i], simplified_points_list[(i + 1) % num_points]) for i
                       in range(num_points)]
            average_length = np.mean(lengths)
            threshold_length = SIMPLIFICATION_THRESHOLD_FACTOR * average_length

            for i in range(num_points):
                # The segment is between point i and point i+1
                if lengths[i] < threshold_length:
                    # Identify the four points involved in the decision
                    p_prev_idx = (i - 1 + num_points) % num_points
                    p0_idx = i
                    p1_idx = (i + 1) % num_points
                    p_next_idx = (i + 2) % num_points

                    p_prev = simplified_points_list[p_prev_idx]
                    p0 = simplified_points_list[p0_idx]
                    p1 = simplified_points_list[p1_idx]
                    p_next = simplified_points_list[p_next_idx]

                    # --- Scenario A: Remove p0 ---
                    # The new neighboring segments would be (p_prev, p1) and (p1, p_next)
                    len_A1 = calculate_distance(p_prev, p1)
                    len_A2 = calculate_distance(p1, p_next)
                    diff_A = abs(len_A1 - len_A2)

                    # --- Scenario B: Remove p1 ---
                    # The new neighboring segments would be (p_prev, p0) and (p0, p_next)
                    len_B1 = calculate_distance(p_prev, p0)
                    len_B2 = calculate_distance(p0, p_next)
                    diff_B = abs(len_B1 - len_B2)

                    # Decide which point to remove
                    if diff_A <= diff_B:
                        # Removing p0 results in more balanced segments
                        del simplified_points_list[p0_idx]
                    else:
                        # Removing p1 results in more balanced segments
                        del simplified_points_list[p1_idx]

                    point_removed_this_pass = True
                    break  # Restart the while loop to re-evaluate with the new point set

            if not point_removed_this_pass:
                break  # No points were removed, so the simplification is done

        # Convert back to a NumPy array for further processing
        final_points = np.array(simplified_points_list)
    else:
        final_points = sorted_points

    # Extract final coordinates for saving and plotting
    if len(final_points) > 0:
        x_coords = final_points[:, 0]
        y_coords = final_points[:, 1]
    else:
        x_coords, y_coords = [], []

    # --- Save the sorted and simplified data to a new CSV file ---
    if len(final_points) > 0:
        sorted_df = pd.DataFrame({
            'x-coordinate': x_coords,
            'y-coordinate': y_coords
        })
        sorted_output_path = 'sorted_interface.csv'
        sorted_df.to_csv(sorted_output_path, index=False)
        print(f"Sorted and simplified coordinates saved to '{sorted_output_path}'")

    # To create a closed loop for plotting, append the first node's coordinates
    if len(final_points) > 1:
        x_line = np.append(x_coords, x_coords[0])
        y_line = np.append(y_coords, y_coords[0])
    else:
        x_line, y_line = x_coords, y_coords

    # --- Plotting ---
    plt.figure(figsize=(10, 8))

    # Plot the lines between nodes (faces)
    plt.plot(x_line, y_line, marker='', linestyle='-', color='royalblue', label='Faces')

    # Plot the dots at the nodes
    plt.plot(x_coords, y_coords, marker='o', linestyle='', color='darkorange', label='Nodes')

    # Add labels and a title for clarity
    plt.xlabel('X-Coordinate')
    plt.ylabel('Y-Coordinate')
    plt.title('Interface Plot (Sorted and Simplified)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)

    # Ensure the aspect ratio is equal to avoid distortion
    plt.axis('equal')

    # Invert y-axis to match typical coordinate systems in some simulations
    # plt.gca().invert_yaxis() # Uncomment if needed

    # Show the plot
    plt.show()

except FileNotFoundError:
    print(f"Error: '{file_path}' not found. Please make sure the file is in the same directory as the script.")
except ValueError as e:
    print(f"Error: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")