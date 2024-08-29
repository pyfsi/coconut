from coconut.tools import flatten_concatenation

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy import integrate

# different cases to be plotted
common_path = '../../pc_development/'
case_paths = ['Stefan_new/case_results.pickle', 'Stefan_new/fine/case_results.pickle', 'Stefan_fixed_itfT/case_results_backup.pickle', 'Stefan_fixed_itfT/fine/case_results.pickle']
legend_entries = ['VarT - coarse', 'VarT - fine', 'CstT - coarse', 'CstT - fine']

itf_faces = [10, 100, 10, 100]
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
line_styles = ['y--', 'g--', 'r--', 'b--', 'k--', 'k-.', 'c--']
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
            L = 0.05  # m, length of original solid domain
            LF = 0.5 - disp_x/L

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
            plt.savefig('figures/itf-disp-x.png')
            plt.show()
            plt.close()
        else:
            # Read Fluent validation files
            read_file_1 = pd.read_csv(r'fluent_val_files/lf_time_coarse.out', delimiter='\s+', skiprows=[0, 1, 2])  # Path of Fluent out-file
            read_file_2 = pd.read_csv(r'fluent_val_files/lf_time.out', delimiter='\s+', skiprows=[0, 1, 2])
            read_file_1.to_csv(r'fluent_val_files/lf_time_coarse.csv', index=None)
            read_file_2.to_csv(r'fluent_val_files/lf_time_fine.csv', index=None)
            data_array_1 = np.loadtxt('fluent_val_files/lf_time_coarse.csv', delimiter=',')
            data_array_2 = np.loadtxt('fluent_val_files/lf_time_fine.csv', delimiter=',')
            LF_val_1 = data_array_1[:, 1]
            time_val_1 = data_array_1[:, 2]
            LF_val_2 = data_array_2[:, 1]
            time_val_2 = data_array_2[:, 2]
            try:
                os.remove("fluent_val_files/lf_time_coarse.csv")
                os.remove("fluent_val_files/lf_time_fine.csv")
            except:
                pass

            time_val_1[0] = 0.0
            time_val_2[0] = 0.0

            # Simple quasi-steady model
            k = 1.5  # W/mK, conduction coefficient of liquid
            dT = 10  # K, the maintained temperature difference over the liquid domain
            rho = 870  # kg/m^3, PCM density
            LH = 179000  # J/kg, latent heat
            B = (k * dT) / (rho * LH)  # m^2/s
            print("B =", B, "m^2/s")

            time_ana = np.linspace(0, 2000, 2001)
            dx = -0.5 * L + 0.5 * np.sqrt(L ** 2 + 4 * B * time_ana)
            LF_ana = 0.5 + dx / L

            line, = plt.plot(time_val_1, LF_val_1, line_styles[len(legend_entries)], label="Fluent - coarse")
            lines_lf.append(line)
            line, = plt.plot(time_val_2, LF_val_2, line_styles[len(legend_entries)+1], label="Fluent - fine")
            lines_lf.append(line)
            line, = plt.plot(time_ana, LF_ana, line_styles[len(legend_entries)+2], label="Ana - steady")
            lines_lf.append(line)

            plt.ylabel('Liquid fraction')
            plt.xlabel('Time [s]')
            plt.legend(handles=lines_lf)
            plt.savefig('figures/liq-frac.png')
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
            line, = plt.plot(time, heat_flux.flatten(), line_styles[j], label=name)
            lines_heat_flux.append(line)

        # Plot interface heat flux in time
        plt.ylabel('Heat flux [W/m^2]')
        plt.xlabel('Time [s]')
        plt.legend(handles=lines_heat_flux)
        plt.savefig('figures/itf-hf.png')
        plt.show()
        plt.close()

# Plot residuals
"""
res = np.array(flatten_concatenation(results['case']['residual']))
it = np.arange(res.size)
line4, = plt.plot(it, res)
plt.ylabel('Residual')
plt.xlabel('Nr. of iterations')
plt.show()
plt.close()
"""