from coconut.coupling_components.component import Component
from coconut.data_structure.variables import variables_dimensions
from coconut import tools

from scipy.spatial import cKDTree
import numpy as np


def create(parameters):
    raise NotImplementedError('this class can only be used as super-class')


class MapperInterpolator(Component):
    def __init__(self, parameters):
        super().__init__()

        # store settings
        self.settings = parameters['settings']
        self.balanced_tree = self.settings.get('balanced_tree', False)
        self.n_nearest = 0  # must be set in sub-class!

        # initialization
        self.n_from, self.n_to = None, None
        self.coords_from, self.coords_to = None, None
        self.tree = None
        self.distances = None
        self.nearest = None
        self.coeffs = None

        # get list with directions
        self.directions = []
        if type(self.settings['directions']) != list:
            raise TypeError('Directions must be a list')
        for direction in self.settings['directions']:
            if direction.lower() not in ['x', 'y', 'z']:
                raise ValueError(f'"{direction}" is not a valid direction')
            if direction.lower() != direction:
                raise ValueError(f'"{direction}" must be lowercase')
            self.directions.append(direction.lower() + '0')
            if len(self.directions) > 3:
                raise ValueError(f'Too many directions given')

        # get scaling (optional)
        self.scaling = None
        if 'scaling' in self.settings:
            self.scaling = np.array(self.settings['scaling']).reshape(1, -1)
            if self.scaling.size != len(self.directions):
                raise ValueError(f'Scaling must have same length as directions')

    def initialize(self, model_part_from, model_part_to):
        super().initialize()

        # get coords_from
        self.n_from = model_part_from.size
        self.coords_from = np.zeros((self.n_from, len(self.directions)))
        for j, direction in enumerate(self.directions):
            self.coords_from[:, j] = getattr(model_part_from, direction)

        # get coords_to
        self.n_to = model_part_to.size
        self.coords_to = np.zeros((self.n_to, len(self.directions)))
        for j, direction in enumerate(self.directions):
            self.coords_to[:, j] = getattr(model_part_to, direction)
        # check if n_from is large enough
        if self.n_from < self.n_nearest:
            tools.print_info(f'Model part {model_part_from} has not enough from-points:'
                             f'{self.n_from} < {self.n_nearest}', layout='warning')
            self.n_nearest = self.n_from

        # check bounding boxes
        tools.check_bounding_box(model_part_from, model_part_to)

        # apply scaling to coordinates
        if self.scaling is not None:
            tools.print_info(f'Scaling {self.scaling} applied for interpolation from {model_part_from.name} '
                             f'to {model_part_to.name}')
            self.coords_from *= self.scaling
            self.coords_to *= self.scaling

        # build and query tree
        if self.balanced_tree:  # time-intensive
            self.tree = cKDTree(self.coords_from)
        else:  # less stable
            self.tree = cKDTree(self.coords_from, balanced_tree=False)

        self.distances, self.nearest = self.tree.query(self.coords_to, k=self.n_nearest)
        self.nearest = self.nearest.reshape(-1, self.n_nearest)

        # check for duplicate points
        self.check_duplicate_points(model_part_from)

    def __call__(self, args_from, args_to):
        # unpack arguments
        interface_from, mp_name_from, var_from = args_from
        interface_to, mp_name_to, var_to = args_to

        # check variables
        if var_from not in variables_dimensions:
            raise NameError(f'Variable "{var_from}" does not exist')
        if var_from != var_to:
            raise TypeError('Variables must be equal')

        # interpolation
        data_from = interface_from.get_variable_data(mp_name_from, var_from)
        data_to = interface_to.get_variable_data(mp_name_to, var_to)
        for i in range(data_to.shape[0]):
            for j in range(data_to.shape[1]):
                data_to[i, j] = np.dot(self.coeffs[i], data_from[self.nearest[i, :], j])
        interface_to.set_variable_data(mp_name_to, var_to, data_to)

    def check_duplicate_points(self, model_part_from):
        # checks only from-points
        tol_warning = 1e-8  # TODO: create optional parameter for this?
        tol_error = 1e-12

        # calculate reference distance (diagonal of bounding box)
        d_ref = np.linalg.norm(self.coords_from.max(axis=0) - self.coords_from.min(axis=0))
        if d_ref == 0.:
            raise ValueError('All from-points coincide')

        # check for duplicate points
        dist, _ = self.tree.query(self.coords_from, k=2)

        msg_1 = f'ModelPart {model_part_from.name}: '
        msg_2 = f' duplicate points found, first duplicate: '

        duplicate = (dist[:, 1] / d_ref < tol_error)
        if duplicate.any():
            raise ValueError(f'{msg_1}{np.sum(duplicate)}{msg_2}' +
                             f'{self.coords_from[np.argmax(duplicate)]}.')

        duplicate = (dist[:, 1] / d_ref < tol_warning)
        if duplicate.any():
            raise Warning(f'{msg_1}{np.sum(duplicate)}{msg_2}' +
                          f'{self.coords_from[np.argmax(duplicate)]}.')
