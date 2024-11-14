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
common_path = '../../pc_development/'
case_paths = ['Stefan_fixed_itfT/case_results.pickle', 'Stefan_fixed_itfT/large_dt/case_results.pickle']
legend_entries = ['Coconut - 0.01 s', 'Coconut - 0.1 s']

itf_faces = [10, 10]
# Plots
disp_plots = False

# load cases
results = {}
for name, path in zip(legend_entries, case_paths):
    with open(os.path.join(common_path, path), 'rb') as file:
        results.update({name: pickle.load(file)})

lines_disp = []
lines_lf = []
lines_heat_flux = []

# make figures
line_styles = ['k--', 'g--', 'r--', 'b--', 'y--', 'k-.', 'c--']
for sol, itf, var, uni in (('solution_x', 'interface_x', 'displacement', 'm'), ('solution_y', 'interface_y', 'heat_flux', 'W/m^2')):
    if var == "displacement":
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
            L = 0.1  # m, length of initial solid domain
            LF = 0.5 + disp_x/L

            if disp_plots:
                line, = plt.plot(time, disp_x.flatten(), line_styles[j], label=name)
                lines_disp.append(line)
            else:
                line, = plt.plot(time, LF.flatten(), line_styles[j], label=name)
                lines_lf.append(line)

        # Plot interface displacement in time
        if disp_plots:
            plt.ylabel('Interface displacement (x) [m]')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_disp)
            plt.savefig('figures/itf-disp-x-Stefan.png')
            plt.show()
            plt.close()
        else:
            # Read Fluent validation files
            read_file_1 = pd.read_csv(r'fluent_val_files/lf_time_coarse.out', delimiter='\s+', skiprows=[0, 1, 2])  # Path of Fluent out-file
            read_file_2 = pd.read_csv(r'fluent_val_files/lf_time_fine.out', delimiter='\s+', skiprows=[0, 1, 2])
            read_file_1.to_csv(r'fluent_val_files/lf_time_coarse.csv', index=None)
            read_file_2.to_csv(r'fluent_val_files/lf_time_fine.csv', index=None)
            data_array_1 = np.loadtxt('fluent_val_files/lf_time_coarse.csv', delimiter=',')
            data_array_2 = np.loadtxt('fluent_val_files/lf_time_fine.csv', delimiter=',')

            LF_val_1 = data_array_1[:, 1]
            time_val_1 = data_array_1[:, 2]
            LF_val_2 = data_array_2[:, 1]
            time_val_2 = data_array_2[:, 2]

            # Arrays for energy balance
            q_in_1 = data_array_1[1:, 3]  # W
            h_avg_1 = data_array_1[1:, 4]  # (J/kg)(kg) = J
            q_in_2 = data_array_2[1:, 3]  # W
            h_avg_2 = data_array_2[1:, 4]  # (J/kg)(kg) = J

            try:
                os.remove("fluent_val_files/lf_time_coarse.csv")
                os.remove("fluent_val_files/lf_time_fine.csv")
            except:
                pass

            time_val_1[0] = 0.0
            time_val_2[0] = 0.0

            # Analytical Stefan solution
            k_l = 1.5  # W/mK, conduction coefficient of liquid
            cp = 2500 # J/kgK, heat capacity
            T_L = 309.15  # K, boundary temperature
            T_m = 299.15  # K, melting temperature
            dT = T_L - T_m  # K, the maintained temperature difference over the liquid domain
            rho = 870  # kg/m^3, PCM density
            LH = 179000  # J/kg, latent heat

            alpha_l = k_l / (rho * cp)
            St_l = (cp * dT) / LH  # Stefan number

            # First: find lambda 'la'
            eq_la = lambda x: x * M.exp(x ** 2) * M.erf(x) - St_l / M.sqrt(M.pi)
            la = fsolve(eq_la, 0.2)[0]

            x_ini = 0.05  # m
            x_end = 0.1  # m, total domain length
            t_ini = x_ini ** 2 / (4 * la ** 2 * alpha_l)
            # t_end = x_end ** 2 / (4 * la ** 2 * alpha_l)
            t_end = t_ini + 600 # solution time

            dt = 0.01
            m = M.ceil((t_end - t_ini)/dt)
            time_ana = np.linspace(t_ini, t_end, m+1)
            x_front_ana = 2 * la * np.sqrt(alpha_l * time_ana)
            LF_ana = x_front_ana/x_end
            q_ana = rho*LH*la*M.sqrt(alpha_l)*time_ana**(-1/2)

            line, = plt.plot(time_val_1, LF_val_1, line_styles[len(legend_entries)], label="Fluent - coarse")
            lines_lf.append(line)
            line, = plt.plot(time_val_2, LF_val_2, line_styles[len(legend_entries)+1], label="Fluent - fine")
            lines_lf.append(line)
            line, = plt.plot(time_ana-t_ini, LF_ana, line_styles[len(legend_entries)+2], label="Ana - stefan")
            lines_lf.append(line)

            plt.ylabel('Liquid fraction')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_lf)
            plt.savefig('figures/liq-frac-Stefan.png')
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
            line, = plt.plot(time[1:-1], -heat_flux.flatten()[1:-1], line_styles[j], label=name)
            lines_heat_flux.append(line)

        line, = plt.plot(time_ana - t_ini, q_ana, line_styles[len(legend_entries)], label="Ana - stefan")
        lines_heat_flux.append(line)

        # Plot interface heat flux in time
        plt.ylabel('Heat flux [W/m^2]')
        plt.xlabel('Time [s]')
        plt.legend(handles=lines_heat_flux)
        plt.savefig('figures/itf-hf-Stefan.png')
        plt.show()
        plt.close()

# Check energy balance of Fluent simulations
if not disp_plots:
    vol = 0.1*0.01*1 # m^3
    sims = ["Coarse", "Fine"]
    q_in = [q_in_1, q_in_2]
    time_val = [time_val_1, time_val_2]
    h_tot = [h_avg_1*rho*vol, h_avg_2*rho*vol]

    for j, sim in enumerate(sims):
        int_flux = integrate.trapezoid(q_in[j], time_val[j][1:])
        delta_h = h_tot[j][-1] - h_tot[j][0]
        diff = int_flux - delta_h
        proc = diff/delta_h*100

        print(sim + ":")
        print('Integrated heat flux is', int_flux, 'J.')
        print('Enthalpy difference is', delta_h, 'J.')
        print('Difference is', diff, 'J, or', proc, '%.')
        print("\n")