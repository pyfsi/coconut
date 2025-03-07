import matplotlib.pyplot as plt
import numpy as np
import math as M
import pandas as pd
from scipy.optimize import fsolve
from scipy.special import erf
from coconut.examples.post_processing.post_processing import *

# Functions
def wall_hf_stefan(k, a, dT, St, t):
    eq_la = lambda la: St / M.sqrt(M.pi) - la * M.exp(la ** 2) * erf(la)
    lam = fsolve(eq_la, 0.2)
    return k * dT / (np.sqrt(M.pi * a * t) * erf(lam))

# Read and load data
sim_file = 'CFD_2/report-file.out'
validation_path = '../2D_post_processing/Alexiades_validation/'

pp = PostProcess('case_results.pickle')
sx = pp.add_subset(interface='interface_x', model_part='boundary_in_nodes')
sy = pp.add_subset(interface='interface_y', model_part='boundary_out_faces')

# Interface location
itf_loc_val = np.loadtxt(validation_path + 'interface_loc.txt', delimiter=' ')
y_val = itf_loc_val[:,0]
x_val_100 = itf_loc_val[:,1]
x_val_200 = itf_loc_val[:,2]
x_val_450 = itf_loc_val[:,3]
x_val_700 = itf_loc_val[:,4]

x_sim_100 = sx.get_values('coordinates', 'x')[10000,:].flatten()
y_sim_100 = sx.get_values('coordinates', 'y')[10000,:].flatten()
x_sim_200 = sx.get_values('coordinates', 'x')[20000,:].flatten()
y_sim_200 = sx.get_values('coordinates', 'y')[20000,:].flatten()
x_sim_450 = sx.get_values('coordinates', 'x')[45000,:].flatten()
y_sim_450 = sx.get_values('coordinates', 'y')[45000,:].flatten()
x_sim_700 = sx.get_values('coordinates', 'x')[70000,:].flatten()
y_sim_700 = sx.get_values('coordinates', 'y')[70000,:].flatten()

line_val_100, = plt.plot(x_val_100, y_val, '*k', label='Validation 100 s')
line_sim_100, = plt.plot(x_sim_100, y_sim_100, '-k', label='Simulation 100 s')
line_val_200, = plt.plot(x_val_200, y_val, '*g', label='Validation 200 s')
line_sim_200, = plt.plot(x_sim_200, y_sim_200, '-g', label='Simulation 200 s')
line_val_450, = plt.plot(x_val_450, y_val, '*r', label='Validation 450 s')
line_sim_450, = plt.plot(x_sim_450, y_sim_450, '-r', label='Simulation 450 s')
line_val_700, = plt.plot(x_val_700, y_val, '*b', label='Validation 700 s')
line_sim_700, = plt.plot(x_sim_700, y_sim_700, '-b', label='Simulation 700 s')
plt.ylabel('y-coordinate [m]')
plt.xlabel('x-coordinate [m]')
plt.xlim((0, 0.1))
plt.title('Validation of explicit simulation: interface location')
plt.legend(handles=[line_val_100, line_sim_100, line_val_200, line_sim_200, line_val_450, line_sim_450, line_val_700, line_sim_700])
plt.savefig('Validation_figures/Validation_itf_loc.png')
plt.show()
plt.close()

# Material properties and problem parameters
k = 60 # W/mK
rho = 7500 # kg/m^3
cp = 200 # J/kgK
alpha = k / (rho * cp) # thermal diffusivity
Tm = 505 # K
T_wall = 508 # K
LH = 60000 # J/kg
Stefan = cp * (T_wall - Tm) / LH # Stefan number
vol_tot = 0.1*0.1*1.0 # m^3
L = 0.1 # m

# Variables at different time instances
var_time_val = np.loadtxt(validation_path + 'variables_in_time.txt', delimiter=' ')
time_val = var_time_val[0,:]
LF_val = var_time_val[1,:]
Nu_val = var_time_val[2,:]*2
stream_max_val = var_time_val[3,:] * 10**(-6)
stream_min_val = var_time_val[4,:] * 10**(-4)
u_max_val = var_time_val[5,:]
u_min_val = var_time_val[6,:]
v_max_val = var_time_val[7,:]
v_min_val = var_time_val[8,:]

# Nusselt plot digitizer
t_pd = [55.87786259541985, 58.62595419847328, 61.37404580152671, 64.12213740458014, 67.78625954198473, 71.45038167938931, 76.03053435114504, 82.44274809160305, 87.02290076335878, 92.51908396946564, 101.6793893129771, 109.00763358778626, 116.33587786259541, 126.41221374045801, 136.4885496183206, 144.73282442748092, 154.8091603053435, 168.54961832061068, 181.3740458015267, 192.36641221374046, 203.3587786259542, 209.7709923664122, 216.18320610687022, 226.2595419847328, 240, 253.74045801526717, 269.3129770992366, 285.80152671755724, 300.45801526717554, 321.5267175572519, 339.8473282442748, 362.7480916030534, 386.5648854961832, 410.381679389313, 433.28244274809157, 460.76335877862596, 469.9236641221374, 475.41984732824426, 479.08396946564886, 481.8320610687023, 490.9923664122137, 499.23664122137404, 515.7251908396946, 533.1297709923664, 554.1984732824427, 582.5954198473282, 599.0839694656488, 616.4885496183206, 628.3969465648855, 642.1374045801526, 656.7938931297709, 674.1984732824427, 689.7709923664122, 699.8473282442748]
Nu_pd = [14.991150442477876, 14.589970501474927, 14.247787610619469, 14, 13.598820058997049, 13.174041297935103, 12.808259587020649, 12.359882005899705, 12.005899705014748, 11.616519174041297, 11.132743362831858, 10.731563421828909, 10.412979351032448, 10, 9.68141592920354, 9.410029498525073, 9.11504424778761, 8.761061946902654, 8.501474926253687, 8.265486725663717, 8.100294985250738, 7.935103244837758, 7.828908554572271, 7.687315634218289, 7.510324483775811, 7.345132743362832, 7.191740412979351, 7.061946902654867, 6.943952802359882, 6.814159292035399, 6.71976401179941, 6.601769911504425, 6.507374631268437, 6.424778761061947, 6.342182890855457, 6.259587020648968, 6.224188790560472, 5.976401179941003, 6.188790560471976, 5.787610619469026, 5.905604719764012, 5.893805309734513, 5.846607669616519, 5.799410029498525, 5.752212389380531, 5.716814159292035, 5.705014749262537, 5.6932153392330385, 5.6932153392330385, 5.6932153392330385, 5.6932153392330385, 5.6932153392330385, 5.705014749262537, 5.705014749262537]
t_pd = np.array(t_pd)
Nu_pd = np.array(Nu_pd)*2

read_file = pd.read_csv(sim_file, delimiter='\s+', skiprows=[0, 1, 2])
read_file.to_csv('Validation_figures/report_file.csv', index=None)
data_sim = np.loadtxt('Validation_figures/report_file.csv', delimiter=',')

time_sim = data_sim[:,10]
LF_sim = data_sim[:,1]/vol_tot
u_max_sim = data_sim[:,2]
u_min_sim = data_sim[:,3]
v_max_sim = data_sim[:,4]
v_min_sim = data_sim[:,5]
stream_max_sim = data_sim[:,6]
stream_min_sim = data_sim[:,7]
q = data_sim[:,8] # Simulation heat flux with convection included
# q_c = wall_hf_stefan(k, alpha, T_wall-Tm, Stefan, time_sim) # Conduction only heat flux according to Stefan solution at the hot wall
# Nu_sim = q/q_c # Definition according to reference -> does not give right results
T_avg = data_sim[:,9]
h = q / (T_wall - T_avg) # heat transfer coefficient
Nu_sim = h*L/k # different Nusselt number definition than reference

# Plots (simulation data starting from 1 s)
# Liquid fraction
line_sim, = plt.plot(time_sim[100:], LF_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, LF_val, '*r', label='Validation')
plt.ylabel('Liquid fraction [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_LF.png')
plt.show()
plt.close()

# Nusselt number
line_sim, = plt.plot(time_sim[100:], Nu_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, Nu_val, '*r', label='Validation')
plt.ylabel('Nusselt number [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_Nu.png')
plt.show()
plt.close()

# Nusselt number from plot digitizer
line_sim, = plt.plot(time_sim[100:], Nu_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(t_pd, Nu_pd, '-*r', label='Val. plot dig.')
plt.ylabel('Nusselt number [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Plot_digitizer_Nu.png')
plt.show()
plt.close()

# Maximum stream function
line_sim, = plt.plot(time_sim[100:], stream_max_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, stream_max_val, '*r', label='Validation')
plt.ylabel('Max stream function [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_stream_max.png')
plt.show()
plt.close()

# Minimum stream function
line_sim, = plt.plot(time_sim[100:], stream_min_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, stream_min_val, '*r', label='Validation')
plt.ylabel('Min stream function [-]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_stream_min.png')
plt.show()
plt.close()

# Maximum x-velocity
line_sim, = plt.plot(time_sim[100:], u_max_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, u_max_val, '*r', label='Validation')
plt.ylabel('Max x-velocity [m/s]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_u_max.png')
plt.show()
plt.close()

# Minimum x-velocity
line_sim, = plt.plot(time_sim[100:], u_min_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, u_min_val, '*r', label='Validation')
plt.ylabel('Min x-velocity [m/s]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_u_min.png')
plt.show()
plt.close()

# Maximum y-velocity
line_sim, = plt.plot(time_sim[100:], v_max_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, v_max_val, '*r', label='Validation')
plt.ylabel('Max y-velocity [m/s]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_v_max.png')
plt.show()
plt.close()

# Minimum y-velocity
line_sim, = plt.plot(time_sim[100:], v_min_sim[100:], '-k', label='Simulation')
line_val, = plt.plot(time_val, v_min_val, '*r', label='Validation')
plt.ylabel('Min y-velocity [m/s]')
plt.xlabel('Time [s]')
plt.legend(handles=[line_val, line_sim])
plt.savefig('Validation_figures/Validation_v_min.png')
plt.show()
plt.close()

# Transverse profiles
transverse_val = np.loadtxt(validation_path + 'transverse_profiles.txt', delimiter=' ')
y_trans_val = transverse_val[:,0]
#stream_val_200 = transverse_val[:,1]
#stream_val_450 = transverse_val[:,3]
#stream_val_700 = transverse_val[:,5]
#T_val_200 = transverse_val[:,2]
#T_val_450 = transverse_val[:,4]
#T_val_700 = transverse_val[:,6]

fluent_profiles = ['stream_200.txt', 'temp_200.txt', 'stream_450.txt', 'temp_450.txt', 'stream_700.txt', 'temp_700.txt']
x_label = ['Stream function at t = 200 s', 'Temperature at t = 200 s', 'Stream function at t = 450 s', 'Temperature at t = 450 s', 'Stream function at t = 700 s', 'Temperature at t = 700 s']
for i, profile_name in enumerate(fluent_profiles):
    read_file = pd.read_csv('CFD_2/' + profile_name, delimiter='\s+', skiprows=[0, 1, 2, 3])
    read_file.to_csv('Validation_figures/' + profile_name.replace('.txt', '.csv'), index=None)
    profile_sim = np.loadtxt('Validation_figures/' + profile_name.replace('.txt', '.csv'), delimiter=',')

    sorted_ind = np.argsort(profile_sim[:,0])

    line_sim, = plt.plot(profile_sim[:,1][sorted_ind], profile_sim[:,0][sorted_ind], '-k', label='Simulation')
    line_val, = plt.plot(transverse_val[:,i+1], y_trans_val, '--*r', label='Validation')
    plt.ylabel('y [m]')
    plt.xlabel(x_label[i])
    plt.legend(handles=[line_val, line_sim])
    plt.savefig('Validation_figures/Val_profile_' + profile_name.replace('.txt', '') + '.png')
    plt.show()
    plt.close()