import matplotlib
import matplotlib.patches as patches
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
matplotlib.rcParams['figure.figsize'] = 20, 8

t0 = time.time()
pp = PostProcess('../case_results.pickle')
pp.print_summary()
print('_____')
pp.print_info()
print('_____')
print(pp)
print('_____')

TS = 5.0e-6
T = [1000, 2000, 3000, 4000]

sx = pp.add_subset(interface='interface_x')
print(sx.get_values('coordinates').shape)

plt.figure()
ax = plt.gca()
x1, x2, y1, y2 = 0.0, 0.3, -0.0015, 0.0015
axins = ax.inset_axes([0.3, 0.67, 0.67, 0.3], xlim=(x1, x2), ylim=(y1, y2), xticks=[0.0, 0.1, 0.2, 0.3],
                      yticks=[-0.0015, 0.0, 0.0015])

# plot nozzle contours
balloon_upper = patches.Arc((0.0048, 0.04064), 0.08, 0.08, angle=0.0, theta1=180, theta2=270, ec='grey', fill=False,
                            lw=2, ls='--')
balloon_lower = patches.Arc((0.0048, -0.04064), 0.08, 0.08, angle=0.0, theta1=90, theta2=180, ec='grey', fill=False,
                            lw=2, ls='--')
ax.add_patch(balloon_upper)
ax.add_patch(balloon_lower)
balloon_upper = patches.Arc((0.0048, 0.04064), 0.08, 0.08, angle=0.0, theta1=180, theta2=270, ec='grey', fill=False,
                            lw=2, ls='--')
balloon_lower = patches.Arc((0.0048, -0.04064), 0.08, 0.08, angle=0.0, theta1=90, theta2=180, ec='grey', fill=False,
                            lw=2, ls='--')
axins.add_patch(balloon_upper)
axins.add_patch(balloon_lower)
plt.plot([0.0048, 0.0352, 0.04045, 0.283102382], [0.00064, 0.00064, 0.00088, 0.0013], color='grey', lw=2, ls='--')
plt.plot([0.0048, 0.0352, 0.04045, 0.283102382], [-0.00064, -0.00064, -0.00088, -0.0013], color='grey', lw=2, ls='--',
         label='Contact body')
axins.plot([0.0048, 0.0352, 0.04045, 0.283102382], [0.00064, 0.00064, 0.00088, 0.0013], color='grey', lw=2, ls='--')
axins.plot([0.0048, 0.0352, 0.04045, 0.283102382], [-0.00064, -0.00064, -0.00088, -0.0013], color='grey', lw=2, ls='--')

A0 = sx.get_values('coordinates')[0]
plt.plot(A0[:, 0], A0[:, 1], label=f't = {0 * TS * 1000.0:.0f} ms', ls=':', color='k')
axins.plot(A0[:, 0], A0[:, 1], ls=':', color='k')

for t in T:
    A0 = sx.get_values('coordinates')[t]
    plt.plot(A0[:, 0], A0[:, 1], label=f't = {t * TS * 1000.0:.0f} ms')
    axins.plot(A0[:, 0], A0[:, 1])

ax.indicate_inset_zoom(axins, edgecolor='black')
plt.xlabel('x-coordinate [m]')
plt.ylabel('y-coordinate [m]')
plt.tight_layout()
plt.subplots_adjust(bottom=0.225)
plt.figlegend(loc='lower center', ncol=6, handles=ax.get_legend_handles_labels()[0])

t1 = time.time()
print(f'Elapsed time: {(t1 - t0):.2f} s')

plt.show()
