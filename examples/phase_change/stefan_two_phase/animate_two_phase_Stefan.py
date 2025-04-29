from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy import integrate
import math as M
from scipy.optimize import fsolve
from scipy.special import erf
from scipy import integrate

# different cases to be plotted
common_path = './'
case_paths = ['Fluent_v2023R1_smoothing/case_results.pickle']
legend_entries = ['new']
fluent_dir = 'fluent_validation/'
fluent_cases = ['lf_fine.out', 'lf_coarse.out']
fluent_legend = ['Fluent - fine', 'Fluent - coarse']

itf_faces = [100]

# Settings
disp_plots = False # plot displacement instead of liquid fraction
fluent_val = False # Is Fluent validation data available?

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

lines_disp = []
lines_lf = []
lines_heat_flux = []
arrays = []
lines_error = []

# make figures
line_styles = ['b--', 'g--', 'r--', 'c--', 'y--', 'k--', 'k-.']
for sol, itf, var, uni in (('solution_x', 'interface_x', 'displacement', 'm'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    if var == "displacement":
        if not disp_plots:
            # Analytical Stefan solution (part 1)
            k_s = 0.024  # W/mK, conduction coefficient of solid
            k_l = 1.5  # W/mK, conduction coefficient of liquid
            cp = 2500  # J/kgK, heat capacity
            T_S = 289.15  # K, boundary temperature solid
            T_L = 309.15  # K, boundary temperature liquid
            T_m = 299.15  # K, melting temperature
            rho = 870  # kg/m^3, PCM density
            LH = 179000  # J/kg, latent heat
            L = 0.1  # m, length of initial solid domain

            alpha_l = k_l / (rho * cp)
            alpha_s = k_s / (rho * cp)
            St_l = cp * (T_L - T_m) / LH # Stefan number in liquid domain
            St_s = cp * (T_m - T_S) / LH # Stefan number in solid domain
            nu = M.sqrt(alpha_l / alpha_s)

            # First: find lambda 'la'
            eq_la = lambda la: St_l / (M.exp(la ** 2) * M.erf(la)) - St_s / (nu * M.exp(nu ** 2 * la ** 2) * M.erfc(nu * la)) - la * M.sqrt(M.pi)
            la = fsolve(eq_la, 0.2)[0]

            x_ini = 0.05  # m
            x_end = L  # m, total domain length
            t_ini = x_ini ** 2 / (4 * la ** 2 * alpha_l)
            # t_end = x_end ** 2 / (4 * la ** 2 * alpha_l)

        for j, name in enumerate(legend_entries):
            solution = results[name][sol]
            interface = results[name][itf]
            dt = results[name]['delta_t']
            time_step_start = results[name]['timestep_start']

            # Store interface data for displacement and LF plot
            disp_x = np.zeros((1, np.shape(solution)[1]))
            for f in [0 + i*3 for i in range((np.shape(solution)[0]-itf_faces[j])//3)]:
                disp_x += solution[f,:]
            disp_x /= (np.shape(solution)[0]-itf_faces[j])//3
            time = time_step_start + dt*np.array(range(np.shape(solution)[1]))

            # Plot liquid fraction
            LF = 0.5 + disp_x/L

            # determine max time reached by Coconut simulations
            if j == 0:
                t_min_co = time[-1]
                t_max_co = time[-1]
                min_index = len(time) - 1
            else:
                if time[-1] < t_min_co:
                    t_min_co = time[-1]
                    min_index = len(time) - 1
                if time[-1] > t_max_co:
                    t_max_co = time[-1]

            if disp_plots:
                line, = plt.plot(time, disp_x.flatten(), line_styles[j], label=name)
                lines_disp.append(line)
            else:
                line, = plt.plot(time + t_ini, LF.flatten(), line_styles[j], label=name)
                lines_lf.append(line)
                arrays.append(LF.flatten())

        # Plot interface displacement in time
        if disp_plots:
            plt.ylabel('Interface displacement (x) [m]')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_disp)
            plt.show()
            plt.close()
        else:
            # Read Fluent validation files
            if fluent_val:
                for j, name in enumerate(fluent_legend):
                    file = common_path + fluent_dir + fluent_cases[j] # Path of Fluent out-file
                    read_file = pd.read_csv(file, delimiter='\s+', skiprows=[0, 1, 2])
                    read_file.to_csv(common_path + fluent_dir + 'report_file.csv', index=None)
                    data_array = np.loadtxt(common_path + fluent_dir + 'report_file.csv', delimiter=',')

                    LF_val = data_array[:, 1]
                    time_val = data_array[:, 2]

                    # Arrays for energy balance
                    q_in = data_array[1:, 3]  # W
                    q_out = data_array[1:, 4]  # W
                    h_avg = data_array[1:, 5]  # (J/kg)(kg) = J

                    try:
                        os.remove(common_path + fluent_dir + 'report_file.csv')
                    except:
                        pass

                    # determine max time reached by Fluent simulations
                    if j == 0:
                        t_max_fl = time_val[-1]
                    else:
                        if time_val[-1] > t_max_fl:
                            t_max_fl = time_val[-1]

                    # set correct initial conditions
                    time_val[0] = 0.0
                    if 'full' in fluent_cases[j]:
                        LF_val[0] = 0.0
                    else:
                        time_val += t_ini

                    # plot fluent results
                    line, = plt.plot(time_val, LF_val, line_styles[len(legend_entries)+j], label=fluent_legend[j])
                    lines_lf.append(line)

                    # Check energy balance
                    vol = 0.1 * 0.01 * 1  # m^3
                    h_tot = h_avg * rho * vol
                    q_tot = q_in + q_out
                    int_flux = integrate.trapezoid(q_tot, time_val[1:])
                    delta_h = h_tot[-1] - h_tot[0]
                    diff = int_flux - delta_h
                    proc = diff / delta_h * 100

                    print(name + ":")
                    print('Integrated heat flux is', int_flux, 'J.')
                    print('Enthalpy difference is', delta_h, 'J.')
                    print('Difference is', diff, 'J, or', proc, '%.')
                    print("\n")
            else:
                t_max_fl = 0

            # Stefan solution (part 2)
            t_end = t_ini + max(t_max_co, t_max_fl)  # solution time
            t_end_error = t_ini + max(t_min_co, t_max_fl)

            dt = 0.01
            m = M.ceil((t_end - t_ini) / dt)
            time_ana = np.linspace(t_ini, t_end, m + 1)
            x_front_ana = 2 * la * np.sqrt(alpha_l * time_ana)
            LF_ana = x_front_ana / x_end
            q_ana = rho * LH * la * M.sqrt(alpha_l) * time_ana ** (-1 / 2)
            start_index = 0
            end_index = M.ceil((t_end_error - t_ini) / dt)

            # plot analytical
            line, = plt.plot(time_ana, LF_ana, line_styles[len(legend_entries) + len(fluent_legend)], label="Ana - stefan")
            lines_lf.append(line)

            # Plot liquid fraction
            plt.ylabel('Liquid fraction')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_lf)
            # plt.tight_layout()
            # plt.savefig('legend_lf.png', dpi=150)
            plt.show()
            plt.close()

            # calculate error of coconut simulations -> only works if all time step sizes are equal!
            for i, lf_sim in enumerate(arrays):
                min_index = max(end_index-start_index, min_index)
                error = np.abs(LF_ana[start_index+1:end_index] - lf_sim[1:min_index])#/np.abs(LF_ana[start_index+1:end_index] - LF_ana[0])
                line, = plt.plot(time_ana[start_index+1:end_index], error, line_styles[i], label=legend_entries[i])
                lines_error.append(line)

            # Log-plot of simulation error
            plt.ylabel('Liquid fraction - error')
            plt.yscale('log')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_error)
            plt.show()
            plt.close()

    if var == "heat_flux":
        for j, name in enumerate(legend_entries):
            solution = results[name][sol]
            interface = results[name][itf]
            dt = results[name]['delta_t']
            time_step_start = results[name]['timestep_start']

            heat_flux = solution[-1,:]
            time = time_step_start + dt * np.array(range(np.shape(solution)[1]))
            line, = plt.plot(t_ini + time[1:-1], -heat_flux.flatten()[1:-1], line_styles[j], label=name)
            lines_heat_flux.append(line)

        line, = plt.plot(time_ana, q_ana, line_styles[len(legend_entries)], label="Ana - stefan")
        lines_heat_flux.append(line)

        # Plot interface heat flux in time
        plt.ylabel('Heat flux [W/m^2]')
        plt.xlabel('Time [s]')
        plt.legend(handles=lines_heat_flux)
        plt.show()
        plt.close()