from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

# different cases to be plotted
common_path = ''
case_paths = ['case_results.pickle']
legend_entries = ['fsi2']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    results.update({name: pickle.load(open(os.path.join(common_path, path), 'rb'))})

# make figure and create animation for each case
animation_figure_displacement = AnimationFigure()  # figure for displacement animations
animation_figure_pressure = AnimationFigure()  # figure for pressure animations
animation_figure_traction_x = AnimationFigure()  # figure for pressure animations
animation_figure_traction_y = AnimationFigure()  # figure for pressure animations
animation_figure_coordinates = AnimationFigure()  # figure for coordinate animations
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni, ani_fig, comp in (("solution_x", "interface_x", "displacement", "m",
                                           animation_figure_displacement, 1),
                                          ("solution_y", "interface_y", "pressure", "Pa",
                                           animation_figure_pressure, 1),
                                          ("solution_y", "interface_y", "traction", "Pa",
                                           animation_figure_traction_x, 0),
                                          ("solution_y", "interface_y", "traction", "Pa",
                                           animation_figure_traction_y, 1),
                                          ("solution_x", "interface_x", "coordinates", "m",
                                           animation_figure_coordinates, 1)):
    for j, name in enumerate(legend_entries):
        for mp in [m['model_part'] for m in results[name][itf].parameters]:
            solution = results[name][sol]
            interface = results[name][itf]
            dt = results[name]["delta_t"]
            time_step_start = results[name]["timestep_start"]
            # create animate object
            animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var, name=name,
                                              model_part_name=mp)
            # select points and component of variable to plot
            coordinates = animation.coordinates

            mask_x = (coordinates[:, 0] > -np.inf)
            mask_y = (coordinates[:, 1] > 0)
            mask_z = (coordinates[:, 2] >= 0)
            abscissa = 0  # x-axis
            component = comp

            animation.initialize(mask_x, mask_y, mask_z, abscissa, component)
            animation.line.set_color(colors[j % len(colors)])
            animation.line.set_linestyle(line_styles[0])
            animation.line.set_marker('o')
            animation.line.set_markersize(2)

    ani_fig.figure.axes[0].set_ylabel(f"{var} ({uni})")
    if var == 'traction':
        ani_fig.figure.axes[0].set_ylabel(f"{'x' if comp == 0 else 'y'}-{var} ({uni})")
    ani_fig.figure.axes[0].set_xlabel("x-coordinate (m)")
    # ani_fig.figure.axes[0].legend()
    ani_fig.figure.tight_layout()
    # or make figure active using plt.figure(ani_fig.number) and use plt.xlabel("") type commands etc.

circle = plt.Circle((0.2, 0.2), 0.05, color='tab:blue', fill=False)
animation_figure_coordinates.figure.axes[0].add_patch(circle)
animation_figure_coordinates.figure.axes[0].set_xlim((0.1, 0.8))
animation_figure_coordinates.figure.axes[0].set_ylim((0.1, 0.3))
animation_figure_coordinates.figure.axes[0].set_aspect('equal', adjustable='box')
animation_figure_coordinates.figure.axes[0].set_title('Turek benchmark FSI2')
animation_figure_coordinates.figure.set_figwidth(6)
animation_figure_coordinates.figure.set_figheight(2.5)
animation_figure_coordinates.figure.set_dpi(300)
animation_figure_coordinates.figure.tight_layout()


for animation_figure, name in ((animation_figure_displacement, 'displacement'),
                               (animation_figure_pressure, 'pressure'),
                               (animation_figure_traction_x, 'traction_x'),
                               (animation_figure_traction_y, 'traction_y'),
                               (animation_figure_coordinates, 'coordinates')):
    animation_figure.make_animation(frames=range(animation_figure.time_step_start,
                                                 animation_figure.time_step_start + animation_figure.time_steps + 1,
                                                 1), #10
                                    interval=1)
    # animation_figure_pressure.make_plot(50)

    save = False
    movie_name = f"{name}.mp4"
    if save:
        # set up formatting for the movie files: mp4-file
        plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
        writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)

        animation_figure.animation.save(movie_name, writer=writer)

plt.show()
