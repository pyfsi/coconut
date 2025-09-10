from odbAccess import *
import os
import sys
import numpy as np

odb = openOdb('CSM_Time2000.odb', readOnly=True)

step = odb.steps.values()[0]

write_dir = 'forces_deformed_geometry'
if not os.path.exists(write_dir):
    os.mkdir(write_dir)

frame = step.frames[1]

yarn = odb.rootAssembly.instances['INSTANCE-YARN']
output_data_nodal = np.zeros((len(yarn.nodes), 13))
displacement_nodal = frame.fieldOutputs['U'].getSubset(region=yarn)
for j, node in enumerate(yarn.nodes):
    instance_nodeID = node.label
    output_nodeID = displacement_nodal.values[j].nodeLabel
    if instance_nodeID != output_nodeID:
        print('Iteration ' + str(j) + ': Instance node label differs from output node label!')
        print('yarn.node.label: ' + str(instance_nodeID) + '; fo_values.nodeLabel: ' + str(output_nodeID))
        sys.stdout.flush()
        sys.exit(1)
    disp = np.array(displacement_nodal.values[j].dataDouble)
    coords = np.array(node.coordinates)
    output_data_nodal[instance_nodeID - 1, 0] = instance_nodeID
    output_data_nodal[instance_nodeID - 1, 1:4] = coords
    output_data_nodal[instance_nodeID - 1, 4:7] = disp

force_nodal = frame.fieldOutputs['SF'].getSubset(region=yarn).getSubset(ELEMENT_NODAL)
strain_nodal = frame.fieldOutputs['SE'].getSubset(region=yarn).getSubset(ELEMENT_NODAL)
for i in range(len(force_nodal.values)):
    tensor = force_nodal.values[i].data
    node = force_nodal.values[i].nodeLabel
    if 1 < node < (len(yarn.elements) + 1):
        # on nodes connecting elements: two values (one for each neighbouring element)
        output_data_nodal[node - 1, 7:10] += tensor / 2.
    else:  # on midpoint node, or first or last yarn node: no averaging
        output_data_nodal[node - 1, 7:10] += tensor

for i in range(len(strain_nodal.values)):
    tensor = strain_nodal.values[i].data
    node = strain_nodal.values[i].nodeLabel
    if 1 < node < (len(yarn.elements) + 1):
        # on nodes connecting elements: two values (one for each neighbouring element)
        output_data_nodal[node - 1, 10:] += tensor / 2.
    else:  # on midpoint node, or first or last yarn node: no averaging
        output_data_nodal[node - 1, 10:] += tensor

np.savetxt(os.path.join(write_dir, 'SF_deformed_yarn_timestep2000_nodes.txt'), output_data_nodal,
           header='nodeID, x0, y0, z0, dx, dy, dz, SF1, SF2, SF3, SE1, SE2, SE3',
           fmt=('%15i', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e', '%15.6e'))
