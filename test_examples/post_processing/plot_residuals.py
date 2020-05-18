import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

# This script plots the convergence residuals of the coupling iterations in every time step up to a cut-off tolerance.
# To disable the cut-off tolerance set it at 0.
# To generate a result file, include a boolean "save_results" in the settings of the coupled solver with value True.
# Give a name to the case by including the string {"name": "a_name"} in the settings of the coupled solver.

tolerance = 1e-8  # cut-off tolerance

# different cases to be plotted
common_path = "../../test_examples/"
case_names = ["results"]
case_paths = ["tube_tube_flow_tube_structure/results"]

# load cases
results = {}
for name, path in zip(case_names, case_paths):
    results.update({name: pickle.load(open(os.path.join(common_path, path), 'rb'))})

# reference case
case_reference = case_names[0]


def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]


def to_tolerance(residuals, tolerance):
    """Cuts off list of residual arrays at one value below tolerance."""
    residuals_out = []
    for ls in residuals:
        i = 0
        ls_out = []
        while i < len(ls) - 1 and ls[i] > tolerance:
            ls_out.append(ls[i])
            i += 1
        ls_out.append(ls[i])
        residuals_out.append(ls_out)
    return residuals_out


# find maximum number of iterations per time step over all cases
it_max = 0
residual_list = []  # list containing residual lists corresponding for the different cases
# each residual list contains a nested list: [ts0[it0, it1, ...], ts1[it0, ...], ...]
for case in case_names:
    residual_list.append(to_tolerance(results[case]["residual"], tolerance))
    for ls in residual_list[-1]:
        if len(ls) > it_max:
            it_max = len(ls)

# make figure
plt.figure()
residual = None
for case, residuals in zip(case_names, residual_list):
    residual = np.zeros((len(residuals), it_max))
    for i, ls in enumerate(residuals):
        residual[i, :len(ls)] = np.array(ls)
    plt.plot(zero_to_nan(residual.flatten()), 'o-', label=case)
plt.plot(tolerance * np.ones_like(residual.flatten()), 'k')
plt.yscale("log")
plt.ylabel("norm of residual")
plt.xlabel("iteration")
plt.legend()

# average number of iteration per time step
for case, residuals in zip(case_names, residual_list):
    iterations = []
    for ls in residuals:
        iterations.append(len(ls))
    avg_iterations = np.array(iterations).mean()
    print(f"{case}: average # iterations = {avg_iterations:0.2f}")
    print('\t', iterations)

plt.show()
