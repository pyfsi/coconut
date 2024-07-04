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
v_01 = (mu/rho)*(0.332*b*(k_f/k_s)*Pr**(1/3))**(-2)*0.01
v_05 = (mu/rho)*(0.332*b*(k_f/k_s)*Pr**(1/3))**(-2)*0.05
v_10 = (mu/rho)*(0.332*b*(k_f/k_s)*Pr**(1/3))**(-2)*0.10
v_15 = (mu/rho)*(0.332*b*(k_f/k_s)*Pr**(1/3))**(-2)*0.15
print("Velocities = " + str(v_01) + " m/s, " + str(v_05) + " m/s, " + str(v_10) + " m/s, " + str(v_15) + " m/s\n")

# Luikov BL solution
v_in = 12 # m/s
x = np.linspace(0.001, 0.2, 200)
Re = rho*v_in*x/mu
Re_L = rho*v_in*L/mu
Bi = 0.332*b*(k_f/k_s)*Pr**(1/3)*Re**(1/2)/x
Bi_L = 0.664*(b/L)*(k_f/k_s)*Pr**(1/3)*Re_L**(1/2)
Bi_avg = Bi_L*np.ones(np.size(x))
print('Velocity = ' + str(v_in) + " m/s:\n")
print('Avg. Bi = ' + str(Bi_L) + "\n")
print('Min. Bi = ' + str(np.min(Bi)) + "\n")
print("Suggested beta TFFB < " + str(1+(np.min(Bi)-1)/(np.min(Bi)+1)) + "\n")
print('Max. Bi = ' + str(np.max(Bi)) + "\n")
print("Suggested beta FFTB < " + str(1-(np.max(Bi)-1)/(np.max(Bi)+1)) + "\n")

# plot
line1, = plt.plot(x, Bi, 'k-', label='Bi')
line2, = plt.plot(x, Bi_avg, 'k--', label='Bi = 1')
plt.ylabel('Bi [-]')
plt.xlabel('Interface location x [m]')
plt.legend()
plt.savefig('luikov_12.svg')
plt.show()
plt.close()