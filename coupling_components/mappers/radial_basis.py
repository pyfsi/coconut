from coconut.coupling_components.mappers.interpolator import MapperInterpolator
from coconut.coupling_components import tools

from scipy.spatial import distance
from scipy.linalg import solve
import numpy as np
from multiprocessing import Pool, cpu_count


def Create(parameters):
    return MapperRadialBasis(parameters)


# Class MapperRadialBasis: radial basis function interpolation.
class MapperRadialBasis(MapperInterpolator):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.coeffs = None
        self.cond = 0

        # check and store settings
        self.parallel = self.settings['parallel'].GetBool() if self.settings.Has('parallel') else False
        self.shape_parameter = self.settings['shape_parameter'].GetInt() if self.settings.Has('shape_parameter') \
            else 200
        if self.shape_parameter < 2:
            tools.Print(f'Shape parameter is {self.shape_parameter} < 2\n', layout='warning')

        # determine number of nearest neighbours
        if len(self.directions) == 3:
            self.n_nearest = 81
        else:
            self.n_nearest = 9

    def Initialize(self, model_part_from, model_part_to):
        super().Initialize(model_part_from, model_part_to)

        # calculate coefficients
        iterable = []
        for i_to in range(self.n_to):
            nearest = self.nearest[i_to, :]
            iterable.append((self.distances[i_to, :],
                             self.coords_from[nearest, :]))

        if self.parallel:
            processes = cpu_count()
            with Pool(processes=processes) as pool:
                # optimal chunksize automatically calculated
                out = pool.starmap(self.get_coeffs, iterable)
            self.coeffs = np.vstack(tuple(out))
        else:
            self.coeffs = np.zeros(self.nearest.shape)
            for i_to, tup in enumerate(iterable):
                self.coeffs[i_to, :] = self.get_coeffs(*tup).flatten()

        # check condition number
        if self.cond > 1e13:
            tools.Print(f'The highest condition number of the interpolation matrices is {self.cond:.2e} > 1e13\n'
                        f'Decrease the shape parameter to decrease the condition number', layout='warning')

    def get_coeffs(self, distances, coords_from):
        def phi(r):
            return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

        d_ref = distances[-1] * self.shape_parameter

        # create column Phi_to, based on distances to from-points
        d_to = distances.reshape(-1, 1)
        Phi_to = phi(d_to / d_ref)

        # create matrix Phi, based on distances between from-points
        d = distance.squareform(distance.pdist(coords_from))
        Phi = phi(d / d_ref)

        # calculate condition number
        cond = np.linalg.cond(Phi)
        if cond > self.cond:
            self.cond = cond

        # solve system Phi^T c = Phi_t for c (Phi is symmetric)
        coeffs = solve(Phi, Phi_to, sym_pos=True)

        return coeffs.reshape(1, -1)
