from coconut.coupling_components.component import Component
from coconut import data_structure
from scipy.spatial import cKDTree

import numpy as np
import copy


def Create(parameters):
    return MapperAxisymmetric2DToWedge3D(parameters)


# Class MapperAxisymmetric3DTo2D: map 3D to 2D axisymmetric.
class MapperAxisymmetric2DToWedge3D(Component):
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

        if not forward:
            # TODO :Adapt code to take the centerline problem into account: print warning if n_from is odd
            self.n_to = model_part_in.NumberOfNodes()
            self.n_from = self.n_to // 2 #2 stands for divide the nodes from 3D to 2D

            self.coords_to = np.zeros((self.n_to, 3))
            for i, node in enumerate(model_part_in.Nodes):
                self.coords_to[i] = np.array([node.X0, node.Y0, node.Z0])
            self.tree = cKDTree(self.coords_to, balanced_tree=False)
            # print(self.tree.data)

                # for j, direction in enumerate(self.directions):
                #     self.coords_to[i, j] = getattr(node, direction)


            model = data_structure.Model()
            model_part_out = model.CreateModelPart('no_name')
            model_part_out._ModelPart__hist_variables = model_part_in._ModelPart__hist_variables

            self.nearest = np.zeros((self.n_to,)).astype(int)

            i_from = 0

            for i_to,node in enumerate(model_part_in.Nodes):
                coords = np.array([node.X0, node.Y0, node.Z0])
                coords_bis=coords.copy()
                r = coords[self.dir_r]
                z = coords[self.dir_wedge]
                if coords[self.dir_wedge] > 0:
                    coords[self.dir_r] = np.cos(self.angle) * r + np.sin(self.angle) * z
                    coords[self.dir_wedge] = 0

                    model_part_out.CreateNewNode(i_from, *tuple(coords))
                    self.nearest[i_to] = i_from

                    coords_bis[self.dir_wedge]=-coords_bis[self.dir_wedge]
                    dd, ii = self.tree.query(coords_bis, k=1)
                    self.nearest[ii] = i_from

                    i_from += 1

            return model_part_out

        else:
            raise NotImplementedError('Backward Initialization not implemented for Mapper2DAxisymmetricToWedge3D.')

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
                node_to.SetSolutionStepValue(var_to, 0, hist_var_from[self.nearest[i_to]])

        # vector interpolation
        elif var_from.Type() == 'Array':
            hist_var_from = np.zeros((self.n_from, 3))
            for i_from, node_from in enumerate(model_part_from.Nodes):
                hist_var_from[i_from] = node_from.GetSolutionStepValue(var_from)

                # print("hist_from")
                # print(hist_var_from[i_from])

            for i_to, node_to in enumerate(model_part_to.Nodes):
                hist_var_to = hist_var_from[self.nearest[i_to]].copy()


                tmp = hist_var_to[self.dir_r]
                # print("tmp")
                # print(tmp)
                # print(self.dir_r)
                hist_var_to[self.dir_r] = tmp*np.cos(self.angle)

                if node_to.Z0 < 0:
                    hist_var_to[self.dir_wedge] = -tmp*np.sin(self.angle)
                else:
                    hist_var_to[self.dir_wedge] = tmp*np.sin(self.angle)
                # print("hist")
                # print(hist_var_to)
                node_to.SetSolutionStepValue(var_to, 0, hist_var_to)



        # other types of Variables
        else:
            raise NotImplementedError(f'Mapping not yet implemented for Variable of Type {var_from.Type()}.')
