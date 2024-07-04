from coconut.examples.post_processing.animate import AnimationFigure

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

debug = True
reverse = False

# different cases to be plotted
if reverse:
    common_path = '../Transient/reverse/'
else:
    common_path = '../Transient/residual_test/'
f = '/case_results.pickle'

if debug:
    case_paths = ['TFFB_relaxation' + f]
    legend_entries = ['TFFB_relaxation']
else:
    case_paths = ['TFFB_relaxation' + f, 'FFTB_relaxation' + f, 'TFFB_aitken' + f, 'FFTB_aitken' + f, 'TFFB_iqni' + f, 'FFTB_iqni' + f]
    legend_entries = ['TFFB_relaxation', 'FFTB_relaxation', 'TFFB_aitken', 'FFTB_aitken', 'TFFB_iqni', 'FFTB_iqni']

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

# make figure and create animation for each case
animation_figure_temperature = AnimationFigure()  # figure for temperature animations
animation_figure_heat_flux = AnimationFigure()  # figure for heat flux animations
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for sol, itf, var, uni, ani_fig in (('solution_x', 'interface_x', 'temperature', 'K', animation_figure_temperature),
                                    ('solution_y', 'interface_y', 'heat_flux', 'W/m^2', animation_figure_heat_flux)):
    for j, name in enumerate(legend_entries):
        solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]['delta_t']
        time_step_start = results[name]['timestep_start']
        # create animate object
        animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var, name=name)
        # select points and component of variable to plot
        coordinates = animation.coordinates

        mask_x = (coordinates[:, 0] > -0.0005)
        mask_y = (coordinates[:, 1] > -0.0005)
        # mask_z = (abs(coordinates[:, 2]) < 1e-16)  # for nodes (displacement)
        # mask_z = (coordinates[:, 2] > 0) & (coordinates[:, 2] < 0.0005)  # for face centers (pressure, traction)
        mask_z = (coordinates[:, 2] > -0.0005) & (coordinates[:, 2] < 0.0005)  # for both
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

animation_figure_temperature.make_animation()
animation_figure_heat_flux.make_animation()
#animation_figure_temperature.make_plot(50)

save = False
animation_figure = animation_figure_temperature
movie_name = './figures/temperature.mp4'
if save:
    # set up formatting for the movie files: mp4-file
    plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
    writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)

    animation_figure.animation.save(movie_name, writer=writer)

plt.show()
