import numpy as np
import matplotlib.pyplot as plt


def calculate_trajectories_with_domain(filename, p1_start, p2_start):
    """
    Calculates trajectories using robust differentiation and Heun's method.
    Plots the trajectories against the specified domain boundaries.
    """

    # --- 1. Load Data ---
    try:
        data = np.loadtxt(filename, comments='#')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # --- 2. Fix Corrupted Data (Extrapolation) ---
    # Discard corrupted first row
    valid_data = data[1:, :]

    # Extrapolate backwards to recreate t=0 state
    row_1 = valid_data[0, :]
    row_2 = valid_data[1, :]

    dt_initial = row_2[0] - row_1[0]

    # Linear Extrapolation: y0 = 2*y1 - y2
    extrapolated_row = 2 * row_1 - row_2
    extrapolated_row[0] = row_1[0] - dt_initial  # Ensure time is correct

    # Stack the new row on top
    data = np.vstack([extrapolated_row, valid_data])
    print(f"Data fixed. Extrapolated start time: t={data[0, 0]:.4f}s")

    # --- 3. Extract columns ---
    t = data[:, 0]
    cg_x = data[:, 1]
    cg_y = data[:, 2]
    v_x = data[:, 3]
    v_y = data[:, 4]
    theta_deg = data[:, 5]

    theta_rad = np.radians(theta_deg)
    omega = np.gradient(theta_rad, t)

    # --- 4. Integration Function (Heun's Method) ---
    def integrate_path_heun(start_x, start_y):
        path_x = np.zeros(len(t))
        path_y = np.zeros(len(t))
        path_x[0] = start_x
        path_y[0] = start_y

        for i in range(len(t) - 1):
            dt = t[i + 1] - t[i]

            # Step A: Current Velocity
            r_x_curr = path_x[i] - cg_x[i]
            r_y_curr = path_y[i] - cg_y[i]

            v_rot_x_curr = -omega[i] * r_y_curr
            v_rot_y_curr = omega[i] * r_x_curr

            vel_x_curr = v_x[i] + v_rot_x_curr
            vel_y_curr = v_y[i] + v_rot_y_curr

            # Step B: Predictor
            pred_x = path_x[i] + vel_x_curr * dt
            pred_y = path_y[i] + vel_y_curr * dt

            # Step C: Future Velocity
            r_x_next = pred_x - cg_x[i + 1]
            r_y_next = pred_y - cg_y[i + 1]

            v_rot_x_next = -omega[i + 1] * r_y_next
            v_rot_y_next = omega[i + 1] * r_x_next

            vel_x_next = v_x[i + 1] + v_rot_x_next
            vel_y_next = v_y[i + 1] + v_rot_y_next

            # Step D: Corrector
            path_x[i + 1] = path_x[i] + (vel_x_curr + vel_x_next) / 2.0 * dt
            path_y[i + 1] = path_y[i] + (vel_y_curr + vel_y_next) / 2.0 * dt

        return path_x, path_y

    # Compute Trajectories
    p1_x, p1_y = integrate_path_heun(p1_start[0], p1_start[1])
    p2_x, p2_y = integrate_path_heun(p2_start[0], p2_start[1])

    # --- 5. Plotting ---
    plt.figure(figsize=(10, 10))

    # A. Plot Domain Boundaries
    # Domain: 120 deg sector (-30 to 90), Radius 51.8mm
    R_domain = 51.8e-3  # Convert mm to meters
    theta_domain = np.linspace(np.radians(-30), np.radians(90), 200)

    # 1. Arc points
    dom_x = R_domain * np.cos(theta_domain)
    dom_y = R_domain * np.sin(theta_domain)

    # 2. Close the loop (Center -> Start -> Arc -> End -> Center)
    domain_x = np.concatenate(([0], dom_x, [0]))
    domain_y = np.concatenate(([0], dom_y, [0]))

    plt.plot(domain_x, domain_y, 'k-', linewidth=2.5, label='Domain Boundary')
    # Optional: Lightly fill the domain to make it visible
    plt.fill(domain_x, domain_y, 'k', alpha=0.05)

    # B. Plot Trajectories
    plt.plot(cg_x, cg_y, 'k--', alpha=0.4, linewidth=1, label='Center of Mass (CG)')

    plt.plot(p1_x, p1_y, 'r-', linewidth=2, label='Point 1')
    plt.plot(p2_x, p2_y, 'b-', linewidth=2, label='Point 2')

    # C. Markers
    plt.scatter([p1_x[0], p2_x[0]], [p1_y[0], p2_y[0]], c='g', s=60, zorder=5, label='Start')
    plt.scatter([p1_x[-1], p2_x[-1]], [p1_y[-1], p2_y[-1]], c='k', marker='x', s=60, zorder=5, label='End')

    # Formatting
    plt.xlabel('X Position (m)')
    plt.ylabel('Y Position (m)')
    plt.title(f'Trajectory within Domain\n{filename.split("/")[-1]}')
    plt.axis('equal')  # Crucial for circular domain
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


# --- Configuration ---
filename = 'Full_7/CFD_2/rigid-body-report-file.out'

# Initial Coordinates
point_1 = (0.03138, 0.02917)
point_2 = (0.03715, 0.01578)

if __name__ == "__main__":
    calculate_trajectories_with_domain(filename, point_1, point_2)