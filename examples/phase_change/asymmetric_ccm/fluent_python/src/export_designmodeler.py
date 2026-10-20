import csv
import os
from pathlib import Path
import numpy as np

# Resolve base simulation directory
CWD = Path.cwd()
if (CWD / "CFD_1" / "interface_nodes.csv").exists():
    BASE_DIR = CWD
else:
    BASE_DIR = Path(__file__).resolve().parent.parent

INPUT_PATH = BASE_DIR / "CFD_1" / "interface_nodes.csv"
OUTPUT_DIR = BASE_DIR / "manual_interventions"
OUTPUT_PATH = OUTPUT_DIR / "interface_nodes_designmodeler.txt"

GROUP_ID = 1
INITIAL_TOLERANCE_FACTOR = 0.001
DESIGNMODELER_NODE_LIMIT = 1999  # Strictly below 2000 nodes for DesignModeler


def ramer_douglas_peucker(points: np.ndarray, epsilon: float) -> np.ndarray:
    """
    Iterative Ramer-Douglas-Peucker line simplification.
    Preserves points exceeding perpendicular distance epsilon from chord.
    """
    n = len(points)
    if n < 3:
        return points

    keep = np.zeros(n, dtype=bool)
    keep[0] = True
    keep[-1] = True

    stack = [(0, n - 1)]

    while stack:
        start_idx, end_idx = stack.pop()
        if end_idx <= start_idx + 1:
            continue

        start = points[start_idx]
        end = points[end_idx]
        segment = points[start_idx + 1 : end_idx]

        line_vec = end - start
        line_len = np.hypot(line_vec[0], line_vec[1])

        if line_len == 0.0:
            dists = np.hypot(segment[:, 0] - start[0], segment[:, 1] - start[1])
        else:
            dists = np.abs(
                line_vec[0] * (start[1] - segment[:, 1]) -
                line_vec[1] * (start[0] - segment[:, 0])
            ) / line_len

        max_idx_local = np.argmax(dists)
        max_dist = dists[max_idx_local]

        if max_dist > epsilon:
            split_idx = start_idx + 1 + max_idx_local
            keep[split_idx] = True
            stack.append((start_idx, split_idx))
            stack.append((split_idx, end_idx))

    return points[keep]


def main():
    print("=== Converting Interface Nodes for ANSYS DesignModeler ===")
    print(f"  Input:  {INPUT_PATH}")
    print(f"  Output: {OUTPUT_PATH}")

    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Read coordinates
    points = []
    with open(INPUT_PATH, mode="r", newline="") as infile:
        csv_reader = csv.reader(infile)
        next(csv_reader)  # Skip header
        for row in csv_reader:
            if len(row) >= 2:
                points.append((float(row[0].strip()), float(row[1].strip())))

    if len(points) < 3:
        raise ValueError(f"Insufficient points in {INPUT_PATH} to generate a curve.")

    # 2. Nearest-neighbour path sorting
    points = np.array(points)
    n_points = len(points)
    path = [0]
    unvisited = np.ones(n_points, dtype=bool)
    unvisited[0] = False

    for _ in range(n_points - 1):
        last_point = points[path[-1]]
        deltas = points[unvisited] - last_point
        dist_sq = np.sum(deltas ** 2, axis=1)

        nearest_local_index = np.argmin(dist_sq)
        unvisited_indices = np.where(unvisited)[0]
        nearest_global_index = unvisited_indices[nearest_local_index]

        path.append(nearest_global_index)
        unvisited[nearest_global_index] = False

    sorted_points = points[path]

    # 3. Calculate baseline spatial scale
    step_vectors = np.diff(sorted_points, axis=0)
    step_lengths = np.hypot(step_vectors[:, 0], step_vectors[:, 1])
    median_spacing = float(np.median(step_lengths))

    # 4. Adaptive RDP Decimation to satisfy DesignModeler limit
    tolerance_factor = INITIAL_TOLERANCE_FACTOR
    epsilon = tolerance_factor * median_spacing
    filtered_points = ramer_douglas_peucker(sorted_points, epsilon)

    adjustment_steps = 0
    while len(filtered_points) > DESIGNMODELER_NODE_LIMIT and len(filtered_points) > 2:
        adjustment_steps += 1
        tolerance_factor *= 1.20
        epsilon = tolerance_factor * median_spacing
        filtered_points = ramer_douglas_peucker(sorted_points, epsilon)

    print(f"\nDecimation Summary:")
    print(f"  Original node count : {len(sorted_points)}")
    print(f"  Final node count    : {len(filtered_points)} (limit: < 2000)")
    print(f"  Final tolerance     : {tolerance_factor:.6f} (scaled x{tolerance_factor / INITIAL_TOLERANCE_FACTOR:.2f} across {adjustment_steps} adjustments)")
    print(f"  Applied epsilon     : {epsilon:.6e}")

    # 5. Export 3D Curve format for DesignModeler
    with open(OUTPUT_PATH, mode="w", newline="") as outfile:
        point_index = 1
        for x, y in filtered_points:
            outfile.write(f"{GROUP_ID} {point_index} {x} {y} 0\n")
            point_index += 1
        outfile.write(f"{GROUP_ID} 0\n")

    print(f"\nSuccessfully generated DesignModeler file at:\n  {OUTPUT_PATH}\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n[ERROR] Failed to convert interface nodes: {e}")
        exit(1)