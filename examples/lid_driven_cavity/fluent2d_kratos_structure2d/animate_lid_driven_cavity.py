from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

# different cases to be plotted
common_path = ''
case_paths = ['case_results.pickle']
legend_entries = ['lid-driven cavity']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

# make figure and create animation for each case
animation_figure_pressure = AnimationFigure()  # figure for pressure animations
animation_figure_coordinates = AnimationFigure()  # figure for coordinate animations
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni, ani_fig in (('solution_y', 'interface_y', 'pressure', 'Pa', animation_figure_pressure),
                                    ('solution_x', 'interface_x', 'coordinates', 'm', animation_figure_coordinates)):
    for j, name in enumerate(legend_entries):
        solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]['delta_t']
        time_step_start = results[name]['timestep_start']
        # create animate object
        animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var, name=name)
        # select points and component of variable to plot
        coordinates = animation.coordinates

        # XY-plane
        mask_x = (coordinates[:, 0] > -np.inf)
        mask_y = (coordinates[:, 1] > -np.inf)
        mask_z = (coordinates[:, 2] > -np.inf)
        abscissa = 0  # x-axis
        component = 1  # y-component

        animation.initialize(mask_x, mask_y, mask_z, abscissa, component)
        animation.line.set_color(colors[j % len(colors)])
        animation.line.set_linestyle(line_styles[0])
        animation.line.set_marker('o')
        animation.line.set_markersize(2)

    ani_fig.figure.axes[0].set_ylabel(f'{var} ({uni})')
    ani_fig.figure.axes[0].set_xlabel('axial coordinate (m)')
    ani_fig.figure.axes[0].legend()
    ani_fig.figure.tight_layout()
    # or make figure active using plt.figure(ani_fig.number) and use plt.xlabel('') type commands etc.

plt.figure(animation_figure_coordinates.number)
plt.plot([0, 0, 1, 1], [0, 1, 1, 0])
plt.ylim([-.1, 1.1])
animation_figure_coordinates.time_position = (0.1, 0.8)
animation_figure_coordinates.figure.axes[0].set_aspect('equal', adjustable='box')
animation_figure_coordinates.figure.set_dpi(300)
animation_figure_coordinates.figure.tight_layout()

frames = None
animation_figure_pressure.make_animation(frames=frames)
animation_figure_coordinates.make_animation(frames=frames)
# animation_figure_pressure.make_plot(50)

save = False
animation_figure = animation_figure_coordinates
movie_name = 'coordinates.mp4'
if save:
    # set up formatting for the movie files: mp4-file
    plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
    writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)

    animation_figure.animation.save(movie_name, writer=writer)

plt.show()
