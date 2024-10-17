from coconut.coupling_components.component import Component
from coconut import tools

from scipy.spatial import cKDTree
import numpy as np


def create(parameters):
    raise NotImplementedError('this class can only be used as super-class')


class MapperInterpolatorBounded(Component):
    def __init__(self, parameters):
        super().__init__()

        # store settings
        self.settings = parameters['settings']
        self.balanced_tree = self.settings.get('balanced_tree', False)
        self.check_bounding_box = self.settings.get('check_bounding_box', True)
        self.n_nearest = 0  # must be set in subclass!

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
            if direction.lower() not in ['x', 'y']:
                raise ValueError(f'"{direction}" is not a valid direction')
            if direction.lower() != direction:
                raise ValueError(f'"{direction}" must be lowercase')
            self.directions.append(direction.lower() + '0')
            if len(self.directions) > 2:
                raise ValueError(f'Too many directions given, this mapper only works in 2D')

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

        print(self.coords_from)
        print('\n')
        print(self.coords_to)

        # check if n_from is large enough
        if self.n_from < self.n_nearest:
            tools.print_info(f'Model part {model_part_from} has not enough from-points:'
                             f'{self.n_from} < {self.n_nearest}', layout='warning')
            self.n_nearest = self.n_from

        # check bounding boxes
        if self.check_bounding_box:
            tools.check_bounding_box(model_part_from, model_part_to, self.directions)

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

    def check_duplicate_points(self, model_part_from):
        # checks only from-points
        tol_warning = 1e-8  # TODO: create optional parameter for this?
        tol_error = 1e-12

        # check is self.coords_from contains more than 1 element
        if self.coords_from.shape[0] > 1:

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
