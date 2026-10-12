from coconut.examples.post_processing.post_processing import *

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

if __name__=="__main__":

    # different cases to be plotted
    common_path = ''
    case_paths = ['case_results.pickle']
    legend_entries = ['CoCoNuT']
    data = {}

    pps = [None] * len(case_paths)
    sxs = [None] * len(case_paths)
    for i, path in enumerate(case_paths):
        cwd = os.getcwd()
        pickel_path = os.path.join(cwd, path)
        pps[i] = pp = PostProcess(pickel_path)
        sxs[i] = sx = pp.add_subset(interface='interface_x')

        # select point to plot
        coordinates = sx.get_all_initial_coordinates()
        sx.select_points(abs(coordinates[:, 0] - 0.0) < 1e-16)

        # displacement of center
        ux = sx.get_values('displacement', 'x')
        uy = sx.get_values('displacement', 'y')
        uz = sx.get_values('displacement', 'z')
        displacement = np.mean((uy*uy + uz*uz)**0.5, axis=1)

        # time array
        time = sx.get_all_times()

        data.update({legend_entries[i]: (time, displacement)})

    # plot
    save = False
    xlim = (0, 0.01)
    ylim = (0.00, 1.3e-4)

    _, ax = plt.subplots(figsize=(10, 7))
    for i, name in enumerate(legend_entries):
        time, displacement = data[name]
        plt.plot(time, displacement, label=name + ' radial displacement',
                linewidth=1.5, marker="o")

    plt.xlabel('Time in s')
    plt.ylabel('Radial displacement in m')
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    ax.tick_params(axis='both', direction='in', pad=8, top=True, right=True)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,1))
    plt.tight_layout()
    plt.legend(loc='upper left')
    if save:
        plt.savefig('comparison_openfoam.png', dpi=300)

    plt.show()
