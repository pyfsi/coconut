import matplotlib
from coconut.examples.post_processing.post_processing import *
import time

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
matplotlib.rcParams['figure.figsize'] = 10, 8

t0 = time.time()
pp = PostProcess('../case_results.pickle')
pp.print_summary()
print('_____')
pp.print_info()
print('_____')
print(pp)
print('_____')

plt.figure()

x_loc = 0.283
TS = 5.0e-6

sx = pp.add_subset(interface='interface_x')
A = sx.get_values('coordinates')
print(A.shape)

v = []
T = []
for t in range(1, sx.get_values('coordinates').shape[0], 20):
    A0 = A[t - 1]
    A1 = A[t]

    ind = np.argmin(np.abs(A0[:, 0] - x_loc))
    if abs(A0[ind, 0] - x_loc) <= 0.01:
        v.append(np.linalg.norm(A1[ind, :] - A0[ind, :]) / TS)
        T.append(t)

T = np.array(T)
plt.plot(T * 5e-03, v, c='k', label='Wool staple-fiber')

plt.xticks(np.arange(0, 26, 5))
plt.xlim(0, 25)
plt.ylim(0, 75)
plt.xlabel('Time [ms]')
plt.ylabel('Velocity [m/s]')
ax = plt.gca()
plt.tight_layout()
plt.subplots_adjust(bottom=0.225)
plt.figlegend(loc='lower center', ncol=2, handles=ax.get_legend_handles_labels()[0])

t1 = time.time()
print(f'Elapsed time: {(t1-t0)/60.:.2f} min')

plt.show()
