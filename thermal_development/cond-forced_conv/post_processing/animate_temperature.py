from coconut.examples.post_processing.animate import AnimationFigure

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

flag = False

# Initial interface temperature
T_ini = 20 # °C

# different cases to be plotted
common_path = '../../../thermal_development/conduction-convection/'
case_paths = ['case_results.pickle']
legend_entries = ['case']

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
for sol, itf, var, uni, ani_fig in (('solution_y', 'interface_y', 'temperature', 'K', animation_figure_temperature),
                                    ('solution_x', 'interface_x', 'heat_flux', 'W/m^2', animation_figure_heat_flux)):
    for j, name in enumerate(legend_entries):
        if var == 'temperature':
            solution = results[name][sol] - 273.15
        else:
            solution = results[name][sol]
        interface = results[name][itf]
        dt = results[name]['delta_t']
        time_step_start = results[name]['timestep_start']
        # create animate object
        # animation = ani_fig.add_animation(solution, interface, dt, time_step_start, variable=var, name=name)
        # select points and component of variable to plot
        # coordinates = animation.coordinates
        """
        mask_x = (coordinates[:, 0] > -np.inf)
        mask_y = (coordinates[:, 1] > 0)
        # mask_z = (abs(coordinates[:, 2]) < 1e-16)  # for nodes (displacement)
        mask_z = (coordinates[:, 2] > -0.0005) & (coordinates[:, 2] < 0.0005)  # for face centers (pressure, traction)
        # mask_z = (coordinates[:, 2] >= 1e-16) & (coordinates[:, 2] < 0.0005)  # for both
        abscissa = 1  # y-axis
        component = 0  # x-component in case solution is vector

        animation.initialize(mask_x, mask_y, mask_z, abscissa, component)
        animation.line.set_color(colors[j % len(colors)])
        animation.line.set_linestyle(line_styles[0])
        animation.line.set_marker('o')
        animation.line.set_markersize(2)
        """
        # Store interface data for temperature plot
        if var == "temperature":
            flag = True
            avg_T = np.array([np.mean(solution[:,i]) for i in range(np.shape(solution)[1])])
            avg_T[0] = T_ini
            time = time_step_start + dt*np.array(range(np.shape(solution)[1]))
    """   
    ani_fig.figure.axes[0].set_ylabel(f'{var} ({uni})')
    ani_fig.figure.axes[0].set_xlabel('height on itf (m)')
    ani_fig.figure.axes[0].legend()
    ani_fig.figure.tight_layout()
    # or make figure active using plt.figure(ani_fig.number) and use plt.xlabel('') type commands etc.
    """
"""
animation_figure_temperature.make_animation()
animation_figure_heat_flux.make_animation()

save = False
animation_figure = animation_figure_temperature
movie_name = 'temperature.mp4'
if save:
    # set up formatting for the movie files: mp4-file
    plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
    writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)

    animation_figure.animation.save(movie_name, writer=writer)

plt.show()
plt.close()
"""
if flag:
    # Plot avg interface temperature in time
    # First, read Fluent validation files
    read_file = pd.read_csv(r'itf-temp-1.out', delimiter='\s+', skiprows=[0, 1, 2]) # Path of Fluent out-file
    read_file.to_csv(r'fluent_validation.csv', index=None)
    data_array = np.loadtxt('fluent_validation.csv', delimiter=',')
    T_val = data_array[:,1] - 273.15
    time_val = data_array[:,2]
    try:
        os.remove("fluent_validation.csv")
    except:
        pass

    # Plot
    line1, = plt.plot(time, avg_T, '*r', label="CoCoNuT")
    line2, = plt.plot(time_val, T_val, '--g', label="Fluent")
    plt.ylabel('Interface temperature [°C]')
    plt.xlabel('Time [s]')
    plt.ylim((np.min(avg_T)-1, np.max(avg_T)+1))
    plt.legend(handles=[line1, line2])
    plt.show()
    plt.close()