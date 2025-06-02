import matplotlib.pyplot as plt
import numpy as np
from coconut.examples.post_processing.post_processing import *
from scipy.optimize import curve_fit

# Define the exponential function to fit
def exponential(t, b):
    return np.exp(b * t)

# different cases to be plotted
common_path = '../'
case_paths = ['Faden_split_3/case_results.pickle', 'Faden_full_17/HPC/case_results.pickle']
legend_entries = ['dont retain', 'retain']
dt = [0.1, 0.1] # s

line_styles = ['r--', 'g--', 'b--', 'k--', 'c--']

prediction = False
lines = []
b_params = []

for j, file in enumerate(case_paths):
    pp = PostProcess(common_path + file)

    sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
    sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

    x = sx.get_values('coordinates', 'x')
    y = sx.get_values('coordinates', 'y')

    # Check if the input arrays are valid
    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise TypeError("Input arrays must be NumPy arrays.")
    if x.shape != y.shape:
        raise ValueError("x and y must have the same shape.")
    if x.ndim != 2:
        raise ValueError("Input arrays must be 2D arrays.")

    N, m = x.shape  # Number of time steps and points

    time = dt[j] * np.arange(N) # Time array
    max_area = np.zeros(N)  # Array to store the maximum distances

    for t in range(N):  # Iterate over each time step
        max_area_t = 0
        for i in range(m-1):
            # Calculate the distance between two neighbouring points
            dist = np.sqrt((x[t, i] - x[t, i+1]) ** 2 + (y[t, i] - y[t, i+1]) ** 2)
            if dist > max_area_t:
                max_area_t = dist
        max_area[t] = max_area_t

    line, = plt.plot(time, max_area, line_styles[j], label=legend_entries[j])
    lines.append(line)

    if prediction:
        # Fit exponential curve through normalised plots
        popt, pcov = curve_fit(exponential, time, max_area/max_area[0], p0=(0.1))  # Initial guess: b=0.1
        b = popt
        b_params.append(b[0])

if prediction:
    # Calculate average of b parameters but exclude outliers
    b_params = np.array(b_params)
    q1 = np.quantile(b_params, 0.25)
    q3 = np.quantile(b_params, 0.75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr

    filtered_b = b_params[(b_params >= lower_bound) & (b_params <= upper_bound)]
    average_b = np.mean(filtered_b)

    full_time = 0.1 * np.linspace(0, 72000)
    line, = plt.plot(full_time, 5e-5 * exponential(full_time, average_b), line_styles[len(case_paths)], label="Exp. prediction")
    lines.append(line)

plt.ylabel('Max. face area [m^2]')
plt.xlabel('Time [s]')
plt.legend(handles=lines)
#plt.savefig('face_area.png')
plt.show()
plt.close()

