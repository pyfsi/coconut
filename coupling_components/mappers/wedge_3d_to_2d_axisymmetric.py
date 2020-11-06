from coconut.coupling_components.component import Component
from coconut import data_structure

import numpy as np
import copy


def Create(parameters):
    return MapperWedge3DToAxisymmetric2D(parameters)


# Class MapperWedge3DToAxisymmetric2D: map wedge 3D to 2D axisymmetric.
class MapperWedge3DToAxisymmetric2D(Component):
    def __init__(self, parameters):
        """
        - should there be both forward and backward initializations? no...
        - should there be a check to see whether input geometry is
            axisymmetrical wrt given directions?
        - take swirl into account? that would really be a
            3D axisymmetrical simulation (only geometry is 2D)...
        """
        super().__init__()

        self.settings = parameters['settings']
        self.interpolator = False

        # get axial and radial directions
        dirs = ['X', 'Y', 'Z']
        self.dir_a = dirs.index(self.settings['direction_axial'].GetString())
        self.dir_r = dirs.index(self.settings['direction_radial'].GetString())
        self.dir_wedge = (set([0, 1 ,2]) - set([self.dir_a, self.dir_r])).pop()
        self.angle=np.radians(2.5)

    def Initialize(self, model_part_in, forward):
        super().Initialize()

        if forward:
            # TODO :Adapt code to take the centerline problem into account: print warning if n_from is odd
            self.n_from = model_part_in.NumberOfNodes()
            self.n_to = self.n_from // 2 #2 stands for divide the nodes from 3D to 2D

            model = data_structure.Model()
            model_part_out = model.CreateModelPart('no_name')
            model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables

            self.nearest = np.zeros((self.n_to,)).astype(int)

            i_to = 0

            for i_from,node in enumerate(model_part_in.Nodes):
                coords = np.array([node.X0, node.Y0, node.Z0])
                r = coords[self.dir_r]
                z = coords[self.dir_wedge]
                if coords[self.dir_wedge] > 0:
                    self.nearest[i_to] = i_from
                    coords[self.dir_r] = np.cos(self.angle) * r + np.sin(self.angle) * z
                    coords[self.dir_wedge] = 0

                    model_part_out.CreateNewNode(i_to, *tuple(coords))
                    i_to += 1

            return model_part_out

        else:
            raise NotImplementedError('Backward Initialization not implemented for MapperWedge3DTo2DAxisymmetric.')

    def __call__(self, args_from, args_to):
        model_part_from, var_from = args_from
        model_part_to, var_to = args_to

        # check if both Variables have same Type
        if var_from.Type() != var_to.Type():
            raise TypeError('Variables to be mapped have different Type.')

        # scalar interpolation
        if var_from.Type() == 'Double':
            hist_var_from = np.zeros(self.n_from)
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)

            for i_to, node_to in enumerate(model_part_to.Nodes):
                hist_var_to=hist_var_from[self.nearest[i_to]]
                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((self.n_from, 3))
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)


            for i_to, node_to in enumerate(model_part_to.Nodes):
                hist_var_to = [0., 0., 0.]
                tmp = hist_var_from[self.nearest[i_to]]
                hist_var_to[self.dir_a] = tmp[self.dir_a]
                hist_var_to[self.dir_r] =tmp[self.dir_r] * np.cos(-self.angle)+tmp[self.dir_wedge] * np.cos(np.pi / 2 - self.angle)

                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)

        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')
