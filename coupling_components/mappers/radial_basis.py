from coconut.coupling_components.mappers.interpolator import MapperInterpolator
from coconut import tools

from scipy.spatial import distance
from scipy.linalg import lstsq, null_space, solve
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
        self.include_polynomial = self.settings.get('include_polynomial', True)

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
            iterable.append((self.coords_to[i_to, :], self.distances[i_to, :], self.coords_from[nearest, :],
                            self.shape_parameter, self.include_polynomial))  # it is faster to include the arguments ..
            # than to make get_coeffs a method of the MapperRadialBasis class

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


def get_coeffs(coords_to, distances, coords_from, shape_parameter, include_polynomial):
    def phi(r):
        return (1 - r) ** 4 * (1 + 4 * r) * np.heaviside(1 - r, 0)

    n = coords_from.shape[0]  # number of from points
    d_ref = distances[-1] * shape_parameter

    # create matrix Phi, based on distances between from-points
    d = distance.squareform(distance.pdist(coords_from))
    Phi = phi(d / d_ref)

    # create column Phi_to, based on distances to from-points
    d_to = distances.reshape(-1, 1)
    Phi_to = phi(d_to / d_ref)

    if include_polynomial:
        # determine scale to improve conditioning number
        scale = np.mean(np.linalg.norm(coords_from, axis=1))

        # create matrix p
        p = np.hstack((np.ones((n, 1)), coords_from / scale))

        # create column p_to
        p_to = np.concatenate((np.array([1]), coords_to / scale)).reshape(-1, 1)

        # determine null space of p
        ns = null_space(p)

        # stack matrices
        if ns.size is not 0:
            # if the from-points are collinear (2d,3d) or coplanar (3d), the linear polynomial is not uniquely defined
            # then, the resulting hyperplane is chosen such that there is no change in value when moving orthogonal
            # to this line/plane
            #
            # if the from-points are collinear (2d,3d) or coplanar (3d), the null space of P is not 0
            # the null space contains the coefficients for the equations such as ax+by+cz+d=0 that represent this
            # line/plane
            # the coefficients a, b, c are the orthogonal directions to this line/plane
            # the requirement is added that the coefficients of the hyperplane are orthogonal to these coefficients
            # i.e., require orthogonality of the normal of the resulting hyperplane (projected on the point space)
            # to directions that are orthogonal to the line/plane formed by the from-points

            ns[0] = 0  # ns contains the directions that are orthogonal to the line/plane formed by the from-points
            ns_dim = ns.shape[1]
            A = np.block([[Phi, p],
                          [p.T, np.zeros((p.shape[1], p.shape[1]))],
                          [np.zeros((ns_dim, n)), ns.T]
                          ])  # not square, full (column) rank
        else:
            A = np.block([[Phi, p], [p.T, np.zeros((p.shape[1], p.shape[1]))]])  # symmetric, full rank
        b = np.block([[Phi_to], [p_to]])
    else:
        A = Phi  # symmetric, full rank
        b = Phi_to

    if A.shape[0] == A.shape[1]:
        # solve system A^T c = b for c (A is symmetric) and truncate
        coeffs = solve(A, b, assume_a='sym')[:n]

        # calculate condition number
        cond = np.linalg.cond(A)
    else:
        # solve system A^T c = b for c, minimal norm solution if A^T is column rank deficit, and truncate
        coeffs, residues, rank, s = lstsq(A.T, b, overwrite_a=True, overwrite_b=True)
        coeffs = coeffs[:n]

        # calculate condition number
        cond = s[0] / s[rank - 1]

    return coeffs.reshape(1, -1), cond
