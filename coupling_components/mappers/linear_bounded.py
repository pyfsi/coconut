from coconut.coupling_components.mappers.interpolator_bounded import MapperInterpolatorBounded
from coconut.data_structure.variables import variables_dimensions

import numpy as np
from scipy.spatial import cKDTree
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
                ll = (self.limits[0], self.limits[1]) # lower left vertex
                ul = (self.limits[0], self.limits[3]) # upper left vertex
                lr = (self.limits[2], self.limits[1]) # lower right vertex
                ur = (self.limits[2], self.limits[3]) # upper right vertex
                self.vertices = [ll, ul, lr, ur]
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

        # Find indices of nodes on domain boundaries
        self.boundary_ind, self.boundary_name = find_boundary_nodes(self.coords_to, self.vertices)
        self.boundary_nodes = self.coords_to[self.boundary_ind]

        if np.shape(self.boundary_nodes)[1] == 2:
            self.boundary_nodes_3D = np.append(self.boundary_nodes, np.zeros((np.shape(self.boundary_nodes)[0], 1)), axis=1)
        elif np.shape(self.boundary_nodes)[1] == 3:
            self.boundary_nodes_3D = self.boundary_nodes

        # build and query tree of nodes
        if self.balanced_tree:  # time-intensive
            self.node_tree = cKDTree(self.coords_to)
        else:  # less stable
            self.node_tree = cKDTree(self.coords_to, balanced_tree=False)

        if np.size(self.boundary_ind) >= 1:
            _, self.neighbour = self.node_tree.query(self.boundary_nodes, k=[2])
            self.neighbour = self.neighbour.reshape(-1, self.n_nearest)
            self.neighbour_nodes = self.coords_to[self.neighbour][0]

            if np.shape(self.neighbour_nodes)[1] == 2:
                self.neighbour_nodes_3D = np.append(self.neighbour_nodes, np.zeros((np.shape(self.neighbour_nodes)[0], 1)), axis=1)
            elif np.shape(self.neighbour_nodes)[1] == 3:
                self.neighbour_nodes_3D = self.neighbour_nodes

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

        # project shifted boundary nodes on domain boundary
        if np.size(self.boundary_ind) >= 1:
            moved_boundary_nodes = self.boundary_nodes_3D + data_to[self.boundary_ind]
            moved_neighbour_nodes = self.neighbour_nodes_3D + data_to[self.neighbour][0]
            for i in range(np.shape(moved_boundary_nodes)[0]):
                params = np.polyfit((moved_boundary_nodes[i, 0], moved_neighbour_nodes[i,0]),
                                        (moved_boundary_nodes[i, 1], moved_neighbour_nodes[i,1]),
                                        deg=1)
                coord = intersect(params, self.boundary_name[i], self.vertices)

                if self.boundary_name[i] == 'top' or self.boundary_name[i] == 'bottom':
                    disp_x = coord[0] - self.boundary_nodes_3D[i, 0]
                    disp_y = 0
                elif self.boundary_name[i] == 'left' or self.boundary_name[i] == 'right':
                    disp_x = 0
                    disp_y = coord[1] - self.boundary_nodes_3D[i, 1]

                data_to[int(self.boundary_ind[i]), 0] = disp_x
                data_to[int(self.boundary_ind[i]), 1] = disp_y

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

def find_boundary_nodes(coords, rect_vertices):
    """
    Finds 2D nodes in an array that lie on the boundary of a rectangular domain.

    Args:
        coords: A NumPy array of shape (n, 3) representing the nodes (z = 0).
        rect_vertices: A list of 4 arrays representing the rectangular domain vertices
                      (lower left, upper left, lower right, upper right).

    Returns:
        A NumPy array containing the indices of the nodes that lie on the boundary.
    """

    # Extract vertex coordinates from the list
    ll_x, ll_y = rect_vertices[0]
    ul_x, ul_y = rect_vertices[1]
    lr_x, lr_y = rect_vertices[2]
    ur_x, ur_y = rect_vertices[3]

    # Check if nodes are within the x and y bounds of the rectangle
    x_within_bounds = (coords[:, 0] >= ll_x) & (coords[:, 0] <= ur_x)
    y_within_bounds = (coords[:, 1] >= ll_y) & (coords[:, 1] <= ul_y)

    # Check if nodes are on the left, right, bottom, or top boundaries
    left_edge = (coords[:, 0] == ll_x) & y_within_bounds
    right_edge = (coords[:, 0] == ur_x) & y_within_bounds
    bottom_edge = (coords[:, 1] == ll_y) & x_within_bounds
    top_edge = (coords[:, 1] == ul_y) & x_within_bounds

    # Combine conditions to find nodes on the boundary
    boundary_indices = np.where(left_edge | right_edge | bottom_edge | top_edge)[0]

    # Determine the corresponding boundary names
    boundary_names = []
    for i in boundary_indices:
        if left_edge[i]:
            boundary_names.append('left')
        elif right_edge[i]:
            boundary_names.append('right')
        elif bottom_edge[i]:
            boundary_names.append('bottom')
        else:
            boundary_names.append('top')

    return boundary_indices, boundary_names

def intersect(param, boundary_name, vertices):
    """
    Finds the intersection point of two lines:
      - One line defined by y = a*x + b
      - Domain boundary line

    Args:
      param: slope and intercept of first line [a b]
      boundary_name: 'left, 'right', 'bottom' or 'top'
      vertices: A list of 4 tuples representing the rectangular domain vertices
                      (lower left, upper left, lower right, upper right).

    Returns:
      The intersection point as a tuple (x, y), or None if there is no intersection.
    """

    if boundary_name == 'left':
        p1, p2 = vertices[0], vertices[1]
    elif boundary_name == 'right':
        p1, p2 = vertices[2], vertices[3]
    elif boundary_name == 'bottom':
        p1, p2 = vertices[0], vertices[2]
    elif boundary_name == 'top':
        p1, p2 = vertices[1], vertices[3]

    param_boundary = np.polyfit((p1[0], p2[0]), (p1[1], p2[1]), deg=1)

    # Check for parallel lines
    if param[0] == param_boundary[0]:
        raise ValueError("Lines are parallel and do not intersect.")

    # Calculate intersection point
    x = (param_boundary[1] - param[1])/(param[0] - param_boundary[0])
    y = param[0] * x + param[1]

    return (x, y)