import numpy as np
import matplotlib.pyplot as plt

# properties
rho = 0.3525 # kg/m^3
cp = 1142.6 # J/kgK
mu = 3.95*10**(-5) # kg/ms
k_f = 0.06808 # W/mK
k_s = 0.2876 # W/mK
Pr = 0.6629

# dimensions
L = 0.2 # m
b = 0.01 # m

# boundary conditions
T_in = 1000 # K
T_b = 600 # K

# velocities to get Bi = 1 @ x = 0.01, 0.05, 0.1 & 0.15
v = [2.4, 12, 24, 36]
lines = []
linestyle = ['b', 'k', 'r', 'g']
x = np.linspace(0.001, 0.2, 200)
Bi_1 = np.ones(np.size(x))
line, = plt.plot(x, Bi_1, 'k--', label='Bi = 1')
lines.append(line)

# Luikov BL solution
for j, v_in in enumerate(v):
    x = np.linspace(0.001, 0.2, 200)
    Re = rho*v_in*x/mu
    Re_L = rho*v_in*L/mu
    Bi = 0.332*b*(k_f/k_s)*Pr**(1/3)*Re**(1/2)/x
    line, = plt.plot(x, Bi, linestyle[j], label='v = ' + str(v_in) + ' m/s')

# plot
plt.ylabel('Bi [-]', fontweight ='bold', fontsize = 12)
plt.xlabel('Interface location x [m]', fontweight ='bold', fontsize = 12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.ylim(0, 3)
plt.savefig('luikov_all_vel.png')
plt.show()
plt.close()