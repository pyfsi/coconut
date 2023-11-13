# basic usage: python forward_append_process_line --output_loc ./ --source_folders drawing1 drawing2 flatPass1 etc.. 

# For a list of cases with the extracted deformation gradients
# provided by the utility extractDeformationGradientSeries
# this creates a list of deformation paths from start to finish in
# a process line.

# The number of deformation paths to be created is based on the starting
# mesh.

import numpy as np
from scipy.spatial import KDTree
from decimal import Decimal

import argparse,glob
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Appends deformation gradients together based on the closest points.')

parser.add_argument('--source_folders', metavar='source_folders', nargs='+',
                    help='The location of the case from which to extract the deformation gradient series',required=True)
parser.add_argument('--output_loc', metavar='output_loc', nargs='+',
                    help='Where to save the data',required=True)

# parse the arguments
args = parser.parse_args()
source_folders = args.source_folders
output_loc = args.output_loc[0]

#print('Linking strain paths in folders:',source_folders)

# read in all the deformation files into memory
# should be okay if we aren't doing extra large cases

full_streamlines = []

start_folder = source_folders[0]
other_passes = source_folders[1:]


# contains the values to be added to each pass to correctly line up the inlets
# the z and y values are set so that the minimum coordinate in the start is set to line up
# with the minimum value
pass_x_offsets = []
pass_y_offsets = []
pass_z_offsets = []

last_pos = 0
prev_y_end_min = 0
prev_z_end_min = 0

all_pass_streamlines = [[] for folder_loc in source_folders]
pass_search_trees = []

for idf,folder_loc in enumerate(source_folders):
    pass_streamlines = glob.glob(folder_loc+'/deformation*')

    # for the average and mins of the starting and ending
    # patches at each pass
    s_start_x = []
    s_start_y = []
    s_start_z = []
    s_end_x = []
    s_end_y = []
    s_end_z = []

    for s in pass_streamlines:
        s_data = np.loadtxt(s,skiprows=1,delimiter=',')

        s_start_x.append(s_data[0,-3])
        s_start_y.append(s_data[0,-2])
        s_start_z.append(s_data[0,-1])

    x_offset = -np.mean(s_start_x) + last_pos
    y_offset = -np.min(s_start_y) + prev_y_end_min
    z_offset = -np.min(s_start_z) + prev_z_end_min

    # now for modifying and appending
    pass_yz_start_points = []
    for s in pass_streamlines:
        s_data = np.loadtxt(s,skiprows=1,delimiter=',')

        s_data[:,-3] += x_offset
#        s_data[:,-2] += y_offset
#        s_data[:,-1] += z_offset

        s_end_x.append(s_data[-1,-3])
        s_end_y.append(s_data[-1,-2])
        s_end_z.append(s_data[-1,-1])

        pass_yz_start_points.append(s_data[0,-2:])

        all_pass_streamlines[idf].append(s_data)

    pass_search_trees.append(KDTree(pass_yz_start_points))
    # pass_x_offsets.append(x_offset)
    # pass_y_offsets.append(y_offset)
    # pass_z_offsets.append(z_offset)

    last_pos = np.mean(s_end_x)
    #prev_y_end_min = np.min(s_end_y)
    #prev_z_end_min = np.min(s_end_z)




# now create the full streamlines
streamlines = [[] for s in all_pass_streamlines[0]]

for ids,s in enumerate(all_pass_streamlines[0]):
    for l in s:
        streamlines[ids].append(l)

for ids,s in enumerate(streamlines):
    for idf,pass_streamlines in enumerate(all_pass_streamlines[1:]):
        query_point = s[-1][-2:]

        dist,index = pass_search_trees[idf+1].query(query_point)

        s2 = pass_streamlines[index]

        for l in s2:
            streamlines[ids].append(l)
        s = s2


        # test to show the stream lines
# for s in streamlines:
#     s_a = np.asarray(s)
#     x = s_a[:,-3]
#     y = s_a[:,-2]
#     z = s_a[:,-1]    

#     plt.plot(x,z)
# plt.show()


# write the streamlines        
for ids,s in enumerate(streamlines):
    with open(output_loc+'/streamline_'+str(ids)+'.dat', 'w') as f:
        f.write('F.xx(),F.xy(),F.xz(),F.yx(),F.yy(),F.yz(),F.zx(),F.zy(),F.zz(),Pos.x(),Pos.y(),Pos.z()\n')
        
        for ids,l in enumerate(s):
            f.write(','.join([str(x) for x in l]) + '\n')




# write the deformation paths for RVEs
for ids,s in enumerate(streamlines):
    with open(output_loc+'/deformationPath_'+str(ids)+'.dat', 'w') as f:
        f.write('(\n')
        
        for ids,l in enumerate(s):
            f.write('    (' + str(ids) + '    (')
            f.write(' '.join([str(x) for x in l[0:9]]) + '))\n')

        f.write(')')

