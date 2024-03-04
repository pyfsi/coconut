from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

# different cases to be plotted
common_path = ''
case_paths = ['case_results.pickle']
legend_entries = ['breaking-dam']

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
for sol, itf, var, uni, ani_fig, abs, comp in (("solution_x", "interface_x", "displacement", "m",
                                                animation_figure_displacement, 1, 0),
                                               ("solution_y", "interface_y", "pressure", "Pa",
                                                animation_figure_pressure, 1, 0),
                                               ("solution_y", "interface_y", "traction", "Pa",
                                                animation_figure_traction_x, 1, 0),
                                               ("solution_y", "interface_y", "traction", "Pa",
                                                animation_figure_traction_y, 1, 1),
                                               ("solution_x", "interface_x", "coordinates", "m",
                                                animation_figure_coordinates, 0, 1)
                                               ):
    for i, name in enumerate(legend_entries):
        for j, mp in enumerate([m['model_part'] for m in results[name][itf].parameters]):
            solution = results[name][sol]
            interface = results[name][itf]
            dt = results[name]["delta_t"]
            time_step_start = results[name]["timestep_start"]
            # create animate object
            animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var,
                                              model_part_name=mp, name=mp)
            mask_x = mask_y = mask_z = np.ones(animation.coordinates[:, 0].size, dtype=bool)
            animation.initialize(mask_x, mask_y, mask_z, abs, comp)
            animation.line.set_color(colors[j % len(colors)])
            animation.line.set_linestyle(line_styles[0])
            animation.line.set_marker('o')
            animation.line.set_markersize(2)

    ani_fig.figure.axes[0].set_ylabel(f"{var} ({uni})")
    if var == 'traction':
        ani_fig.figure.axes[0].set_ylabel(f"{'x' if comp == 0 else 'y'}-{var} ({uni})")
    if var == 'coordinates':
        A = 0.1
        H = 0.14
        L = 0.079
        S = 0.005
        G = 0.0025
        points = np.array(
            [[0, 0], [A, 0], [A, 1.1 * H], [0, 1.1 * H], [0, G + L], [np.nan, np.nan], [-S, G + L], [-3 * A, G + L],
             [-3 * A, 0], [0, 0]])
        ani_fig.figure.axes[0].plot(points[:, 0], points[:, 1], color='black', linewidth=1)
    ani_fig.figure.axes[0].set_xlabel("y-coordinate (m)")
    ani_fig.figure.axes[0].legend()
    ani_fig.figure.tight_layout()
    # or make figure active using plt.figure(ani_fig.number) and use plt.xlabel("") type commands etc.

animation_figure_coordinates.figure.axes[0].set_ylim((-0.1 * H, 1.2 * H))
animation_figure_coordinates.figure.axes[0].set_aspect('equal', adjustable='box')
animation_figure_coordinates.figure.axes[0].set_title('Breaking dam')
animation_figure_coordinates.figure.set_dpi(300)
animation_figure_coordinates.figure.tight_layout()


for animation_figure, name in ((animation_figure_displacement, 'displacement'),
                               (animation_figure_pressure, 'pressure'),
                               (animation_figure_traction_x, 'traction_x'),
                               (animation_figure_traction_y, 'traction_y'),
                               (animation_figure_coordinates, 'coordinates')
                               ):
    animation_figure.make_animation(frames=range(animation_figure.time_step_start,
                                                 animation_figure.time_step_start + animation_figure.time_steps + 1,
                                                 1),
                                    interval=1)
    # animation_figure.make_plot(0)

    save = False
    movie_name = f"{name}.mp4"
    if save:
        # set up formatting for the movie files: mp4-file
        plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
        writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)

        animation_figure.animation.save(movie_name, writer=writer)

plt.show()
