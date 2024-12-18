import os

from coconut.examples.post_processing.post_processing import *

# different cases to be plotted
common_path = '../../tube/'
case_paths = ['fluent3d_abaqus3dcse', 'fluent3d_abaqus3d']
common_file_name = 'case_results.pickle'
legend_entries = ['CSE', 'Old']

# create PostProcess instances
pps = []
for path in case_paths:
    pps.append(PostProcess(os.path.join(common_path, path, common_file_name)))

# create SubSets for each case
sxs = []
sys = []
for pp in pps:
    sxs.append(pp.add_subset(interface='interface_x'))
    sys.append(pp.add_subset(interface='interface_y'))

# select points
for ss in sxs + sys:
    initial_coordinates = ss.get_all_initial_coordinates()
    mask1 = abs(initial_coordinates[:, 2]) < 0.0005
    mask2 = initial_coordinates[:, 1] > 0
    ss.select_points(mask1 & mask2)

# make animations
animations = [Animation2dDisplacement(sxs, 'x', 'y', name=legend_entries),
              # Animation2dCoordinates(sxs, 'x', 'y', name=legend_entries),
              Animation2dPressure(sys, 'x', name=legend_entries),
              Animation2dTraction(sys, 'x', 'x', name=legend_entries),
              ]

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

# save animations
save = False
if save:cp
    # save the animations as mp4; make sure FFmpeg is available (UGent cluster: ml FFmpeg
    # writer = ani.FFMpegFileWriter(metadata=dict(artist='CoCoNuT'), fps=24, bitrate=2000)
    # for animation in animations:
    #     animation.set_writer(writer)
    #     animation.save(f'{animation.figure_name}.mp4')

    # save the animations as gif
    for animation in animations:
        animation.save()

plt.show()
