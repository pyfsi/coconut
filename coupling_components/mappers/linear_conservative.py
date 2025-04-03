from coconut.coupling_components.component import Component
from coconut import tools
from coconut.data_structure.variables import variables_dimensions

import numpy as np
import math
from scipy.spatial import cKDTree
from multiprocessing import Pool, cpu_count

def create(parameters):
    return MapperLinearConservative(parameters)

# Class MapperLinearConservative: linear interpolation with projection of nodes at boundaries + volume conservation at cell level.
class MapperLinearConservative(Component):
    def __init__(self, parameters):
        """
        geometrical calculations:
            P_a = point a
            v_ab = vector from a to b

            subscripts:
                0 = to-point
                1 = 1st from-point
                2 = 2nd from-point
                p = to-point projected on line
                n = normal unit vector
                t = tangential unit vector
        """
        super().__init__()

        # store settings
        self.settings = parameters['settings']
        self.balanced_tree = self.settings.get('balanced_tree', False)
        self.parallel = self.settings.get('parallel', False)
        self.check_bounding_box = self.settings.get('check_bounding_box', True)

        # volume-conserving iterations
        self.tolerance = self.settings.get('volume_tolerance', 1e-14)
        self.max_iterations = self.settings.get('mapper_iterations', 20)
        self.relax = self.settings.get('mapper_relaxation', 0.4)
        self.conservative = self.settings.get('mapping_conservative', True)
        self.projection_order = self.settings.get('projection_order', 1) # 0: no projection, 1: 1st order upwind projection, 2: 2nd order upwind
        self.C = None

        # initialization
        self.n_from, self.n_to = None, None
        self.coords_from, self.coords_to = None, None
        self.tree = None
        self.nearest_faces = None
        self.nearest_nodes = None
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

        # create domain
        if self.projection_order != 0:
            self.domain_type = self.settings.get('domain', None)
            self.limits = self.settings.get('limits', None)

            if self.domain_type == 'rectangular':
                if len(self.limits) == 4:
                    ll = (self.limits[0], self.limits[1])  # lower left vertex
                    ul = (self.limits[0], self.limits[3])  # upper left vertex
                    lr = (self.limits[2], self.limits[1])  # lower right vertex
                    ur = (self.limits[2], self.limits[3])  # upper right vertex
                    self.vertices = [ll, ul, lr, ur]
                else:
                    raise ValueError(
                        'The 2D coordinates of the lower left vertex (0) and upper right vertex (1) are required in the following order: [x0, y0, x1, y1].')
            elif self.domain_type is None:
                raise ValueError(
                    'This mapper requires domain specifications. Choose between following domain types: rectangular.')
            else:
                raise ValueError(
                    f'Given domain type {self.domain_type} is not supported. Choose between following domain types: rectangular.')

        # determine number of nearest
        if len(self.directions) == 2:
            self.n_nearest = 2
        else:
            raise ValueError('This mapper only works in 2D')

    def initialize(self, model_part_from, model_part_to):
        super().initialize()

        # get coords_from
        self.n_from = model_part_from.size
        self.coords_from = np.zeros((self.n_from, len(self.directions)))
        for j, direction in enumerate(self.directions):
            self.coords_from[:, j] = getattr(model_part_from, direction)

        if not self.conservative:
            self.C = np.ones((self.n_from, 1))

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
        if self.check_bounding_box:
            tools.check_bounding_box(model_part_from, model_part_to, self.directions)

        # apply scaling to coordinates
        if self.scaling is not None:
            tools.print_info(f'Scaling {self.scaling} applied for interpolation from {model_part_from.name} '
                             f'to {model_part_to.name}')
            self.coords_from *= self.scaling
            self.coords_to *= self.scaling

        # build and query FROM (faces) tree
        if self.balanced_tree:  # time-intensive
            self.tree = cKDTree(self.coords_from)
        else:  # less stable
            self.tree = cKDTree(self.coords_from, balanced_tree=False)

        _, self.nearest_faces = self.tree.query(self.coords_to, k=self.n_nearest)
        self.nearest_faces = self.nearest_faces.reshape(-1, self.n_nearest)

        # check for duplicate points
        self.check_duplicate_points(model_part_from)

        # build and query TO (nodes) tree
        if self.balanced_tree:  # time-intensive
            self.node_tree = cKDTree(self.coords_to)
        else:  # less stable
            self.node_tree = cKDTree(self.coords_to, balanced_tree=False)

        _, self.nearest_nodes = self.node_tree.query(self.coords_from, k=self.n_nearest)
        self.nearest_nodes = self.nearest_nodes.reshape(-1, self.n_nearest)

        if self.projection_order != 0:
            # Find indices of nodes on domain boundaries
            self.boundary_ind, self.boundary_name = find_boundary_nodes(self.coords_to, self.vertices)
            self.boundary_nodes = self.coords_to[self.boundary_ind]
            self.n_endpoints = np.size(self.boundary_ind)

            if np.size(self.boundary_ind) >= 1:
                _, self.neighbour = self.node_tree.query(self.boundary_nodes, k=[2]) # closest neighbour
                self.neighbour = self.neighbour.reshape(-1, self.n_endpoints)
                self.neighbour_nodes = self.coords_to[self.neighbour][0]
                if self.projection_order == 2:
                    _, self.neighbour_2 = self.node_tree.query(self.boundary_nodes, k=[3]) # 2nd closest neighbour
                    self.neighbour_2 = self.neighbour_2.reshape(-1, self.n_endpoints)
                    self.neighbour_2_nodes = self.coords_to[self.neighbour_2][0]


    def __call__(self, args_from, args_to):
        # unpack arguments
        interface_from, mp_name_from, var_from = args_from
        interface_to, mp_name_to, var_to = args_to

        # check variables
        if var_from not in variables_dimensions:
            raise NameError(f'Variable "{var_from}" does not exist')
        if var_to not in variables_dimensions:
            raise NameError(f'Variable "{var_to}" does not exist')

        # find face & node displacement & position of previous time step
        self.prev_disp_from = interface_from.get_variable_data(mp_name_from, 'prev_disp')[:, 0:len(self.directions)]
        self.prev_coord_from = self.coords_from + self.prev_disp_from
        self.prev_disp_to = interface_to.get_variable_data(mp_name_to, 'prev_disp')[:, 0:len(self.directions)]
        self.prev_coord_to = self.coords_to + self.prev_disp_to

        # update boundary & neighbour nodes
        if self.projection_order != 0:
            self.boundary_nodes = self.prev_coord_to[self.boundary_ind]
            self.neighbour_nodes = self.prev_coord_to[self.neighbour][0]
            if self.projection_order == 2:
                self.neighbour_2_nodes = self.prev_coord_to[self.neighbour_2][0]

        # recalculate coefficients based on new face position
        iterable = []
        for i_to in range(self.n_to):
            nearest = self.nearest_faces[i_to, :]
            iterable.append((self.prev_coord_from[nearest, :], self.prev_coord_to[i_to, :]))

        if self.parallel:
            processes = cpu_count()
            with Pool(processes=processes) as pool:
                # optimal chunksize automatically calculated
                out = pool.starmap(get_coeffs, iterable)
            self.coeffs = np.vstack(tuple(out))
        else:
            self.coeffs = np.zeros(self.nearest_faces.shape)
            for i_to, tup in enumerate(iterable):
                self.coeffs[i_to, :] = get_coeffs(*tup).flatten()

        # get single time step displacements
        data_from = interface_from.get_variable_data(mp_name_from, var_from)[:, 0:len(self.directions)]
        data_to = interface_to.get_variable_data(mp_name_to, var_to)[:, 0:len(self.directions)]

        # initialise multiplication coefficients for face displacements
        if self.conservative is True:
            if self.C is None:
                C = np.ones((self.n_from, 1))
            else:
                C = self.C

            it = 0
            C_prev = np.zeros((self.n_from, 1))

            # calculate cell-based volume
            area = interface_from.get_variable_data(mp_name_from, 'area')
            dx = np.linalg.norm(data_from, axis=1)
            dx = dx.reshape(self.n_from, 1)
            cell_vol = area * dx
            node_vol = np.zeros(np.shape(cell_vol))

            if not np.any(dx): # if all face displacements are zero
                C = np.ones((self.n_from, 1))
            else:
                # perform loop to determine multiplication factor C to compensate for lost or falsely created volume
                while np.linalg.norm(cell_vol - node_vol) > self.tolerance and it < self.max_iterations:
                    it += 1
                    C_prev = C

                    # interpolate and project shifted boundary nodes on domain boundary
                    data_itp = self.interpolate(C * data_from, data_to)

                    # compare cell-based volume to node-based volume
                    for i, c_vol in enumerate(cell_vol):
                        nearest = self.nearest_nodes[i, :]
                        prev_coord = self.prev_coord_to[nearest, :]
                        new_coord = prev_coord + data_itp[nearest, :]
                        x = np.array([prev_coord[0][0], prev_coord[1][0], new_coord[0][0], new_coord[1][0]])
                        y = np.array([prev_coord[0][1], prev_coord[1][1], new_coord[0][1], new_coord[1][1]])
                        node_vol[i] = PolyArea(x, y)

                    C = cell_vol/node_vol
                    C[~np.isfinite(C)] = 1.0
                    C = (1 - self.relax) * C + self.relax * C_prev

                # print warning if C did not converge fast enough
                if it >= self.max_iterations and np.linalg.norm(cell_vol - node_vol) > self.tolerance:
                    warning_text = f'Conservative mapper - C did not converge below tolerance: {np.linalg.norm(cell_vol - node_vol)} > {self.tolerance}.'
                    tools.print_info(warning_text, layout='warning')

            # store converged C
            self.C = C

        # interpolate with final C values
        data_itp = self.interpolate(self.C * data_from, data_to)

        # add z-coordinate to data_itp
        if np.shape(data_itp)[1] == 2:
            data_itp = np.append(data_itp, np.zeros((np.shape(data_itp)[0], 1)), axis=1)

        # Set mapped variable data
        interface_to.set_variable_data(mp_name_to, var_to, data_itp)

    def interpolate(self, data_from, data_to):

        # interpolation
        for i in range(data_to.shape[0]):
            for j in range(data_to.shape[1]):
                data_to[i, j] = np.dot(self.coeffs[i], data_from[self.nearest_faces[i, :], j])

        # boundary projection
        if self.projection_order != 0 and np.size(self.boundary_ind) >= 1:
            moved_boundary_nodes = self.boundary_nodes + data_to[self.boundary_ind]
            moved_neighbour_nodes = self.neighbour_nodes + data_to[self.neighbour][0]
            if self.projection_order == 2:
                moved_neighbour_2_nodes = self.neighbour_2_nodes + data_to[self.neighbour_2][0]
            for i in range(np.shape(moved_boundary_nodes)[0]):
                # list [(x0, x1), (y0, y1)] of 2D faces closest to domain boundaries
                line_coords = [(moved_boundary_nodes[i, 0], moved_neighbour_nodes[i, 0]),
                                    (moved_boundary_nodes[i, 1], moved_neighbour_nodes[i, 1])]
                coord = intersect(line_coords, self.boundary_name[i], self.vertices)

                if self.boundary_name[i] == 'top' or self.boundary_name[i] == 'bottom':
                    disp_x = coord[0] - self.boundary_nodes[i, 0]
                    disp_y = 0
                elif self.boundary_name[i] == 'left' or self.boundary_name[i] == 'right':
                    disp_x = 0
                    disp_y = coord[1] - self.boundary_nodes[i, 1]

                if self.projection_order == 2:
                    line_coords = [(moved_boundary_nodes[i, 0], moved_neighbour_2_nodes[i, 0]),
                                   (moved_boundary_nodes[i, 1], moved_neighbour_2_nodes[i, 1])]
                    coord = intersect(line_coords, self.boundary_name[i], self.vertices)

                    if self.boundary_name[i] == 'top' or self.boundary_name[i] == 'bottom':
                        disp_x_2 = coord[0] - self.boundary_nodes[i, 0]
                        disp_y_2 = 0
                    elif self.boundary_name[i] == 'left' or self.boundary_name[i] == 'right':
                        disp_x_2 = 0
                        disp_y_2 = coord[1] - self.boundary_nodes[i, 1]

                    # 2nd order upwind projection
                    data_to[int(self.boundary_ind[i]), 0] = 1.5 * disp_x - 0.5 * disp_x_2
                    data_to[int(self.boundary_ind[i]), 1] = 1.5 * disp_y - 0.5 * disp_y_2
                else:
                    # 1st order upwind projection
                    data_to[int(self.boundary_ind[i]), 0] = disp_x
                    data_to[int(self.boundary_ind[i]), 1] = disp_y

        return data_to

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

# Functions outside the Interpolator class
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

def intersect(line_coords, boundary_name, vertices):
    """
    Finds the intersection point of two lines:
      - Line between the node at the domain boundary and its nearest node neighbour
      - Domain boundary line

    Args:
      line_coords: list of tuples [(x0, x1), (y0, y1)] with (x0, y0) the node coordinates at the domain boundary and (x1, y1) the neighbouring node coordinates
      boundary_name: 'left, 'right', 'bottom' or 'top'
      vertices: A list of 4 tuples representing the rectangular domain vertices
                      (lower left, upper left, lower right, upper right).

    Returns:
      The intersection point as a tuple (x, y), or None if there is no intersection.
    """

    x_nodes = line_coords[0]
    y_nodes = line_coords[1]
    param, vertical = linear_fit(x_nodes, y_nodes)

    if boundary_name == 'left':
        p1, p2 = vertices[0], vertices[1]
    elif boundary_name == 'right':
        p1, p2 = vertices[2], vertices[3]
    elif boundary_name == 'bottom':
        p1, p2 = vertices[0], vertices[2]
    elif boundary_name == 'top':
        p1, p2 = vertices[1], vertices[3]

    param_boundary, vert_boundary = linear_fit((p1[0], p2[0]), (p1[1], p2[1]))

    if vertical is True and vert_boundary is True:
        raise ValueError("Two vertical lines do not intersect.")
    elif vertical is True and vert_boundary is False:
        x = x_nodes[0]
        y = param_boundary[0] * x + param_boundary[1]
    elif vertical is False and vert_boundary is True:
        x = p1[0]
        y = param[0] * x + param[1]
    else:
        # Check for parallel lines
        if param[0] == param_boundary[0]:
            raise ValueError("Lines are parallel and do not intersect.")

        # Calculate intersection point
        x = (param_boundary[1] - param[1]) / (param[0] - param_boundary[0])
        y = param[0] * x + param[1]

    return (x, y)

def PolyArea(x, y):
    # sort points counter-clockwise
    theta = np.arctan2(y - np.mean(y), x - np.mean(x))
    sorted_index = np.argsort(theta)
    x_sorted = x[sorted_index]
    y_sorted = y[sorted_index]

    # calculate area
    return 0.5 * np.abs(np.dot(x_sorted, np.roll(y_sorted, -1)) - np.dot(y_sorted, np.roll(x_sorted, -1)))

def linear_fit(x, y):
    """
    Finds the coefficents a and b of the line defined by y = a*x + b going through two points

    Args:
      x: tuple with x-coordinates of both points (x1, x2)
      y: tuple with y-coordinates of both points (y1, y2)

    Returns:
      coeff: tuple with coefficients (a, b) such that y = a*x + b goes trough points (x1, y1) and (x2, y2)
      vertical: boolean to indicate whether x1 == x2
    """

    vert = False
    coeff = [0, 0]

    if x[0] != x[1]:
        coeff[0] = (y[1] - y[0]) / (x[1] - x[0])
        coeff[1] = y[0] - coeff[0] * x[0]
    else:
        vert = True

    return np.array(coeff), vert