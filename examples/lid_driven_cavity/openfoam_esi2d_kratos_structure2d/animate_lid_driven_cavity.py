from coconut.examples.post_processing.post_processing import *

# case to be plotted
case_path = 'case_results.pickle'
legend_entry = 'lid-driven cavity'

pp = PostProcess(case_path)

fa_set = dict(skip=9)

# make animations
animations = [Animation2dCoordinates(pp, 'x', 'y', text_location=(0.1, 0.8), func_animation_settings=fa_set),
              Animation2dPressure(pp, 'x', func_animation_settings=fa_set)]

# change figure layout
for i, animation in enumerate(animations):
    # add units to the labels
    ax = animation.get_ax()
    ax.set_xlabel(f'{ax.xaxis.get_label().get_text()} (m)')
    y_label = ax.yaxis.get_label().get_text()
    ax.set_ylabel(f'{y_label} ({"Pa" if "pressure" in y_label or "traction" in y_label else "m"})')

animation_coordinates = animations[0]
animation_coordinates.make_active()
plt.plot([0, 0, 1, 1], [0, 1, 1, 0])
plt.ylim([-.1, 1.1])
animation_coordinates.get_figure().set_dpi(300)
animation_coordinates.get_figure().tight_layout()

# save animations
save = False
if save:
    for animation in animations:
        animation.save()

plt.show()
