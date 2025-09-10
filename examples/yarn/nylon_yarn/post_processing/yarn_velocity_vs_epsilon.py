import matplotlib
from coconut.examples.post_processing.post_processing import *
import time
from os.path import join

matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['lines.linewidth'] = 2
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rcParams['axes.labelpad'] = 15
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['xtick.minor.size'] = 2.4
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rcParams['ytick.minor.size'] = 2.4
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = 20, 8


def get_yarn_velocity(case_name):
    pp = PostProcess(case_name)
    sx = pp.add_subset(interface='interface_x')
    A = sx.get_values('coordinates')
    v = []
    T = []
    for t in range(1, A.shape[0], 20):
        A0 = A[t - 1]
        A1 = A[t]

        ind = np.argmin(np.abs(A0[:, 0] - x_loc))
        if abs(A0[ind, 0] - x_loc) <= 0.01:
            v.append(np.linalg.norm(A1[ind, :] - A0[ind, :]) / TS)
            T.append(t)
    return np.array(v), np.array(T)


t0 = time.time()
plt.figure()
ax = plt.gca()

x_loc = 0.283
TS = 5.0e-6

# Overset case
A = np.loadtxt('yarn_velocity_overset.dat')
plt.plot(A[:, 1], A[:, 2], label=f'Overset', c='k')

# ALM cases
epsilons = ['4R', '3R']
for eps in epsilons:
    velocity, timesteps = get_yarn_velocity(f'results_eps{eps}.pickle')
    plt.plot(timesteps * TS, velocity, label=f'ALM' + r'$\ - \ \epsilon = $' + eps)

plt.xlabel('Time [s]')
plt.ylabel('Velocity [m/s]')
plt.tight_layout()
plt.subplots_adjust(bottom=0.22)
plt.figlegend(loc='lower center', ncol=len(ax.get_legend_handles_labels()[0]), handles=ax.get_legend_handles_labels()[0])

t1 = time.time()
print(f'Elapsed time: {(t1 - t0) / 60.:.2f} min')

plt.show()
