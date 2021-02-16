from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

# different cases to be plotted
common_path = "../../examples/"
case_paths = ["tube_tube_flow_tube_structure/results.pickle"]
legend_entries = ["results"]

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    results.update({name: pickle.load(open(os.path.join(common_path, path), 'rb'))})

# make figure and create animation for each case
animation_figure_displacement = AnimationFigure()  # figure for displacement animations
animation_figure_pressure = AnimationFigure()  # figure for pressure animations
animation_figure_coordinates = AnimationFigure()  # figure for coordinate animations
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni, ani_fig in (("solution_x", "interface_x", "displacement", "m", animation_figure_displacement),
                                    ("solution_y", "interface_y", "pressure", "Pa", animation_figure_pressure),
                                    ("solution_x", "interface_x", "coordinates", "m", animation_figure_coordinates)):
    for j, name in enumerate(legend_entries):
        solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]["delta_t"]
        time_step_start = results[name]["timestep_start"]
        # create animate object
        animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var, name=name)
        # select points and component of variable to plot
        coordinates = animation.coordinates

        # example 1: python solver (YZ-plane)
        python_solver = True
        if python_solver:
            mask_x = (coordinates[:, 0] > -np.inf)
            mask_y = (abs(coordinates[:, 1]) > 0)
            mask_z = (coordinates[:, 2] > -np.inf)
            abscissa = 2  # z-axis
            component = 1  # y-component

        # example 2: fluent solver (XY-plane)
        fluent = False
        if fluent:
            mask_x = (coordinates[:, 0] > -np.inf)
            mask_y = (coordinates[:, 1] > 0)
            mask_z = (abs(coordinates[:, 2]) < 1e-16)
            # mask_z = (coordinates[:, 2] > 0) & (coordinates[:, 2] < 0.0005)
            abscissa = 0  # x-axis
            component = 1  # y-component

        animation.initialize(mask_x, mask_y, mask_z, abscissa, component)
        animation.line.set_color(colors[j % len(colors)])
        animation.line.set_linestyle(line_styles[0])
        animation.line.set_marker('o')
        animation.line.set_markersize(2)

    ani_fig.figure.axes[0].set_ylabel(f"{var} ({uni})")
    ani_fig.figure.axes[0].set_xlabel("axial coordinate (m)")
    ani_fig.figure.axes[0].legend()
    ani_fig.figure.tight_layout()
    # or make figure active using plt.figure(ani_fig.number) and use plt.xlabel("") type commands etc.

animation_figure_displacement.make_animation()
animation_figure_pressure.make_animation()
animation_figure_coordinates.make_animation()
# animation_figure_pressure.make_plot(50)

save = False
animation_figure = animation_figure_displacement
movie_name = "displacement.mp4"
if save:
    # set up formatting for the movie files: mp4-file
    plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
    writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='NicolasDelaissé'), fps=24, bitrate=2000)

    animation_figure.animation.save(movie_name, writer=writer)

plt.show()
