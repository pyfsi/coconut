from coconut.coupling_components.mappers.interpolator import MapperInterpolator
from coconut import tools

from scipy.spatial import distance
from scipy.linalg import solve
import numpy as np
from multiprocessing import Pool, cpu_count


def create(parameters):
    return MapperRadialBasis(parameters)


# Class MapperRadialBasis: radial basis function interpolation.
class MapperRadialBasis(MapperInterpolator):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.coeffs = None

        # check and store settings
        self.parallel = self.settings.get('parallel', False)
        self.shape_parameter = self.settings.get('shape_parameter', 200)
        if self.shape_parameter < 2:
            tools.print_info(f'Shape parameter is {self.shape_parameter} < 2\n', layout='warning')

        # determine number of nearest neighbours
        n_nearest = 81 if len(self.directions) == 3 else 9
        self.n_nearest = self.settings.get('n_nearest', n_nearest)

    def initialize(self, model_part_from, model_part_to):
        super().initialize(model_part_from, model_part_to)

        # calculate coefficients
        iterable = []
        cond = []
        for i_to in range(self.n_to):
            nearest = self.nearest[i_to, :]
            iterable.append((self.distances[i_to, :],
                             self.coords_from[nearest, :],
                             self.shape_parameter))

        if self.parallel:
            processes = cpu_count()
            with Pool(processes=processes) as pool:
                # optimal chunksize automatically calculated
                out = pool.starmap(get_coeffs, iterable)
            self.coeffs = np.vstack(tuple(zip(*out))[0])
            cond = list(zip(*out))[1]
        else:
            self.coeffs = np.zeros(self.nearest.shape)
            for i_to, tup in enumerate(iterable):
                out = get_coeffs(*tup)
                self.coeffs[i_to, :] = out[0].flatten()
                cond.append(out[1])

        # check condition number
        cond = max(cond)
        if cond > 1e13:
            tools.print_info(f'The highest condition number of the interpolation matrices is {cond:.2e} > 1e13\n'
                             f'Decrease the shape parameter (current value is {self.shape_parameter})'
                             f' to decrease the condition number', layout='warning')


def get_coeffs(distances, coords_from, shape_parameter):
    def phi(r):
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

    d_ref = distances[-1] * shape_parameter

    # create column Phi_to, based on distances to from-points
    d_to = distances.reshape(-1, 1)
    Phi_to = phi(d_to / d_ref)

    # create matrix Phi, based on distances between from-points
    d = distance.squareform(distance.pdist(coords_from))
    Phi = phi(d / d_ref)

    # calculate condition number
    cond = np.linalg.cond(Phi)

    # solve system Phi^T c = Phi_t for c (Phi is symmetric)
    coeffs = solve(Phi, Phi_to, sym_pos=True)

    return coeffs.reshape(1, -1), cond
