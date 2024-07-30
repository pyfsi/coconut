import os

from coconut.examples.post_processing.post_processing import *

# different cases to be plotted
common_path = '../../examples/'
case_paths = ['tube/tube_flow_tube_structure/case_results.pickle']
legend_entries = ['tube case']

# create PostProcess instances
pps = []
for path in case_paths:
    pps.append(PostProcess(os.path.join(common_path, path)))

# create SubSets for each case
sxs = []
sys = []
for pp in pps:
    sxs.append(pp.add_subset(interface='interface_x'))
    sys.append(pp.add_subset(interface='interface_y'))

# select points
# not necessary for this case since we want to show all points

# make animations
animations = [Animation2dDisplacement(sxs, 'z', 'y', name=legend_entries),
              Animation2dCoordinates(sxs, 'z', 'y', name=legend_entries),
              Animation2dPressure(sys, 'z', name=legend_entries)]

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
    ax.set_ylabel(f'{y_label} ({"Pa" if "pressure" in y_label else "m"})')

# save animations
save = False
if save:
    # save the animations as mp4; make sure FFmpeg is available (UGent cluster: ml FFmpeg
    # writer = ani.FFMpegFileWriter(metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)
    # for animation in animations:
    #     animation.set_writer(writer)
    #     animation.save(f'{animation.figure_name}.mp4')

    # save the animations as gif
    for animation in animations:
        animation.save()

plt.show()
