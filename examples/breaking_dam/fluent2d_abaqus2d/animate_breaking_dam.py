from coconut.examples.post_processing.post_processing import *

# case to be plotted
case_path = 'case_results.pickle'
legend_entry = 'breaking dam'

pp = PostProcess(case_path)

# make animations
animations = [Animation2dDisplacement(pp, 'y', 'x'),
              Animation2dCoordinates(pp, 'x', 'y'),
              Animation2dPressure(pp, 'y'),
              Animation2dTraction(pp, 'y', 'x', figure_name='Animation2dTractionX'),
              Animation2dTraction(pp, 'y', 'y', figure_name='Animation2dTractionY')]

# change figure layout
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for i, animation in enumerate(animations):
    for j, line in enumerate(animation.get_lines()):
        line.set(color=colors[j % len(colors)], linestyle=line_styles[0], marker='o', markersize=2)
    ax = animation.get_ax()
    ax.legend()  # update legend
    # add units to the labels
    ax.set_xlabel(f'{ax.xaxis.get_label().get_text()} (m)')
    y_label = ax.yaxis.get_label().get_text()
    ax.set_ylabel(f'{y_label} ({"Pa" if "pressure" in y_label or "traction" in y_label else "m"})')

animation_coordinates = animations[1]
A = 0.1
H = 0.14
L = 0.079
S = 0.005
G = 0.0025
points = np.array(
    [[0, 0], [A, 0], [A, 1.1 * H], [0, 1.1 * H], [0, G + L], [np.nan, np.nan], [-S, G + L], [-3 * A, G + L],
     [-3 * A, 0], [0, 0]])
ax = animation_coordinates.get_ax()
ax.plot(points[:, 0], points[:, 1], color='black', linewidth=1)
ax.set_ylim((-0.1 * H, 1.2 * H))
ax.set_title('Breaking dam')

fig = animation_coordinates.get_figure()
fig.set_dpi(300)
fig.tight_layout()

# save animations
save = False
if save:
    for animation in animations:
        animation.save()

plt.show()
