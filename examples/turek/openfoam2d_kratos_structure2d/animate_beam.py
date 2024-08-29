from coconut.examples.post_processing.post_processing import *

# case to be plotted
case_path = 'case_results.pickle'
legend_entry = 'fsi2'

pp = PostProcess(case_path)

fa_set = dict(skip=9)  # func_animation_setting

# make animations
animations = [Animation2dDisplacement(pp, 'x', 'y', func_animation_settings=fa_set),
              Animation2dCoordinates(pp, 'x', 'y', func_animation_settings=fa_set),
              Animation2dPressure(pp, 'x', func_animation_settings=fa_set),
              Animation2dTraction(pp, 'x', 'x', figure_name='Animation2dTractionX', func_animation_settings=fa_set),
              Animation2dTraction(pp, 'x', 'y', figure_name='Animation2dTractionY', func_animation_settings=fa_set)]

# change figure layout
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
line_styles = ['-', '--', ':', '-.']
for i, animation in enumerate(animations):
    for j, line in enumerate(animation.get_lines()):
        line.set(color=colors[j % len(colors)], linestyle=line_styles[j % len(colors)], marker='o', markersize=2)
    ax = animation.get_ax()
    ax.legend()  # update legend
    # add units to the labels
    ax.set_xlabel(f'{ax.xaxis.get_label().get_text()} (m)')
    y_label = ax.yaxis.get_label().get_text()
    ax.set_ylabel(f'{y_label} ({"Pa" if "pressure" in y_label or "traction" in y_label else "m"})')

animation_coordinates = animations[1]
ax = animation_coordinates.get_ax()

circle = plt.Circle((0.2, 0.2), 0.05, color='tab:blue', fill=False)
ax.add_patch(circle)
ax.set_xlim((0.1, 0.8))
ax.set_ylim((0.1, 0.3))
ax.set_title(f'Turek benchmark {legend_entry.capitalize()}')
fig = animation_coordinates.get_figure()
fig.set_figwidth(6)
fig.set_figheight(2.5)
fig.set_dpi(300)
fig.tight_layout()

# save animations
save = False
if save:
    for animation in animations:
        animation.save()

plt.show()
