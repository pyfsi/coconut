from coconut.coupling_components.mappers.interpolator_bounded import MapperInterpolatorBounded
from coconut.data_structure.variables import variables_dimensions

import numpy as np
from multiprocessing import Pool, cpu_count


def create(parameters):
    return MapperLinearBounded(parameters)


# Class MapperLinear: linear interpolation.
class MapperLinearBounded(MapperInterpolatorBounded):
    def __init__(self, parameters):
        """
        geometrical calculations:
            P_a = point a
            v_ab = vector from a to b

            subscripts:
                0 = to-point
                1 = 1st from-point
                2 = 2nd from-point
                3 = 3rd from-point
                p = to-point projected on line/triangle
                n = normal unit vector
                t = tangential unit vector
        """
        super().__init__(parameters)

        # store settings
        self.parallel = self.settings['parallel'] if 'parallel' in self.settings else False
        self.domain_type = self.settings.get('domain', None)
        self.limits = self.settings.get('limits', None)

        # create domain
        if self.domain_type == 'rectangular':
            if len(self.limits) == 4:
                self.ll = np.array(self.limits[0:1]) # lower left vertex
                self.ul = np.array([self.limits[0], self.limits[3]]) # upper left vertex
                self.lr = np.array([self.limits[2], self.limits[1]]) # lower right vertex
                self.ur = np.array(self.limits[2:3]) # upper right vertex
            else:
                raise ValueError('The 2D coordinates of the lower left vertex (0) and upper right vertex (1) are required in the following order: [x0, y0, x1, y1].')
        elif self.domain_type is None:
            raise ValueError('This mapper requires domain specifications. Choose between following domain types: rectangular.')
        else:
            raise ValueError(f'Given domain type {self.domain_type} is not supported. Choose between following domain types: rectangular.')

        # determine number of nearest
        if len(self.directions) == 2:
            self.n_nearest = 2
        else:
            raise ValueError('This mapper only works in 2D')

    def initialize(self, model_part_from, model_part_to):
        super().initialize(model_part_from, model_part_to)

        print(self.coords_from)
        print('\n')
        print(self.coords_to)

    def __call__(self, args_from, args_to):
        # unpack arguments
        interface_from, mp_name_from, var_from = args_from
        interface_to, mp_name_to, var_to = args_to

        # check variables
        if var_from not in variables_dimensions:
            raise NameError(f'Variable "{var_from}" does not exist')
        if var_from != var_to:
            raise TypeError('Variables must be equal')

        # calculate coefficients
        iterable = []
        for i_to in range(self.n_to):
            nearest = self.nearest[i_to, :]
            iterable.append((self.coords_from[nearest, :], self.coords_to[i_to, :]))

        if self.parallel:
            processes = cpu_count()
            with Pool(processes=processes) as pool:
                # optimal chunksize automatically calculated
                out = pool.starmap(get_coeffs, iterable)
            self.coeffs = np.vstack(tuple(out))
        else:
            self.coeffs = np.zeros(self.nearest.shape)
            for i_to, tup in enumerate(iterable):
                self.coeffs[i_to, :] = get_coeffs(*tup).flatten()

        # interpolation
        data_from = interface_from.get_variable_data(mp_name_from, var_from)
        data_to = interface_to.get_variable_data(mp_name_to, var_to)
        for i in range(data_to.shape[0]):
            for j in range(data_to.shape[1]):
                data_to[i, j] = np.dot(self.coeffs[i], data_from[self.nearest[i, :], j])

        # Set mapped varaible data
        interface_to.set_variable_data(mp_name_to, var_to, data_to)

def get_coeffs(coords_from, coord_to):
    coeffs = np.zeros(2)
    P_0 = coord_to
    P_1 = coords_from[0, :]
    P_2 = coords_from[1, :]

    coeffs[0] = line_interpolation_coeff(P_0, P_1, P_2)
    coeffs[1] = 1. - coeffs[0]
    return coeffs.reshape(1, -1)

def line_interpolation_coeff(P_0, P_1, P_2):
    # project P_0 on line
    # *** only necessary if 2D??
    if P_0.size == 1:
        P_p = P_0
    else:
        v_01 = P_1 - P_0
        v_t = (P_2 - P_1) / np.linalg.norm(P_2 - P_1)
        P_p = P_0 + (v_01 - v_t * np.dot(v_01, v_t))

    # check if point lies on line
    if np.dot(P_1 - P_p, P_2 - P_p) <= 0:
        # linear interpolation
        c = np.linalg.norm(P_2 - P_p) / np.linalg.norm(P_2 - P_1)
    else:
        # nearest neighbour interpolation
        c = 1.
    return c