from coconut import data_structure
from coconut.coupling_components.tools import create_instance

import unittest
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def split(coords):
    x, y, z = np.hsplit(coords, 3)
    return x.flatten(), y.flatten(), z.flatten()


class TestMapperInterpolator(unittest.TestCase):
    def setUp(self):
        self.parameters = {'type': 'mappers.nearest',
                           'settings': {'directions': ['z']}}

    # TODO: check tree? check __call__ method?

    def test_instantiation(self):
        self.parameters['settings']['directions'] = ['z']
        mapper = create_instance(self.parameters)
        self.assertListEqual(mapper.directions, ['z0'])

        self.parameters['settings']['directions'] = ['z', 'y', 'x']
        mapper = create_instance(self.parameters)
        self.assertListEqual(mapper.directions, ['z0', 'y0', 'x0'])

        # TODO: add this again, after all capitals have been removed
        # self.parameters['settings'].SetArray('directions', ['Z'])
        # self.assertRaises(ValueError, create_instance, self.parameters)

        self.parameters['settings']['directions'] = ['z0']
        self.assertRaises(ValueError, create_instance, self.parameters)

        self.parameters['settings']['directions'] = ['z', 'y', 'x', 'z']
        self.assertRaises(ValueError, create_instance, self.parameters)

        self.parameters['settings']['directions'] = 'z'
        self.assertRaises(TypeError, create_instance, self.parameters)

    def test_bounding_box_1d(self):
        self.parameters['settings']['directions'] = ['z']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()

        coords = np.array([[0, 0, 0], [0, 0, 1]])
        mp_from = model.create_model_part('mp_from', *split(coords), np.arange(2))
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, 1.01]])))
        mp_to = model.create_model_part('mp_to_2', *split(coords), np.arange(3))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, -.01]])))
        mp_to = model.create_model_part('mp_to_3', *split(coords), np.arange(4))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, 1.1]])))
        mp_to = model.create_model_part('mp_to_4', *split(coords), np.arange(5))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, 1.25]])))
        mp_to = model.create_model_part('mp_to_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, -.25]])))
        mp_to = model.create_model_part('mp_to_6', *split(coords), np.arange(7))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[0, 0, 2], [0, 0, -1]])))
        mp_to = model.create_model_part('mp_to_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

    def test_bounding_box_2d(self):
        self.parameters['settings']['directions'] = ['z', 'x']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()

        coords = np.array([[0, 0, 0], [1, 0, 1]])
        mp_from = model.create_model_part('mp_from', *split(coords), np.arange(2))
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.01, 0, 1.01]])))
        mp_to = model.create_model_part('mp_to_2', *split(coords), np.arange(3))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[-.01, 0, -.01]])))
        mp_to = model.create_model_part('mp_to_3', *split(coords), np.arange(4))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.1, 0, 1.1]])))
        mp_to = model.create_model_part('mp_to_4', *split(coords), np.arange(5))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.25, 0, 1.25]])))
        mp_to = model.create_model_part('mp_to_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[-.25, 0, -.25]])))
        mp_to = model.create_model_part('mp_to_6', *split(coords), np.arange(7))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[2, 0, 2], [-1, 0, -1]])))
        mp_to = model.create_model_part('mp_to_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

    def test_bounding_box_3d(self):
        self.parameters['settings']['directions'] = ['z', 'x', 'y']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()

        coords = np.array([[0, 0, 0], [1, 1, 1]])
        mp_from = model.create_model_part('mp_from', *split(coords), np.arange(2))
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.01, 1.01, 1.01]])))
        mp_to = model.create_model_part('mp_to_2', *split(coords), np.arange(3))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[-.01, -.01, -.01]])))
        mp_to = model.create_model_part('mp_to_3', *split(coords), np.arange(4))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.1, 1.1, 1.1]])))
        mp_to = model.create_model_part('mp_to_4', *split(coords), np.arange(5))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1.25, 1.25, 1.25]])))
        mp_to = model.create_model_part('mp_to_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[-.25, -.25, -.25]])))
        mp_to = model.create_model_part('mp_to_6', *split(coords), np.arange(7))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[2, 2, 2], [-1, -1, -1]])))
        mp_to = model.create_model_part('mp_to_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

    def test_bounding_box_lines(self):
        # check if method works for lines aligned with coordinate axes in 2D
        self.parameters['settings']['directions'] = ['x', 'y']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()
        coords = np.array([[0, 1, 0], [1, 1, 0]])
        mp_from = model.create_model_part('mp_from', *split(coords), np.arange(2))

        coords = np.array([[0, 1, 0], [1.01, 1.01, 0.]])
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.array([[0, 1.05, 0], [1, 1.05, 0]])
        mp_to = model.create_model_part('mp_to_2', *split(coords), np.arange(2))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.array([[0, 1.25, 0], [1, 1.25, 0]])
        mp_to = model.create_model_part('mp_to_3', *split(coords), np.arange(2))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.array([[0, .85, 0], [1, 1.15, 0]])
        mp_to = model.create_model_part('mp_to_4', *split(coords), np.arange(2))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

    def test_bounding_box_planes(self):
        # check if method works for planes aligned with coordinate axes in 3D
        self.parameters['settings']['directions'] = ['x', 'y', 'z']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()
        coords = np.array([[0, 1, 0], [1, 1, 0]])
        mp_from = model.create_model_part('mp_from', *split(coords), np.arange(2))

        coords = np.array([[0, 1, 0], [1.01, 1.01, 0.]])
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))
        mapper.initialize(mp_from, mp_to)
        mapper.finalize()

        coords = np.array([[0, 1.05, 0], [1, 1.05, 0]])
        mp_to = model.create_model_part('mp_to_2', *split(coords), np.arange(2))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.array([[0, 1.25, 0], [1, 1.25, 0]])
        mp_to = model.create_model_part('mp_to_3', *split(coords), np.arange(2))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.array([[0, .85, 0], [1, 1.15, 0]])
        mp_to = model.create_model_part('mp_to_4', *split(coords), np.arange(2))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

    def test_check_duplicate_points(self):
        self.parameters['settings']['directions'] = ['x', 'y', 'z']
        mapper = create_instance(self.parameters)

        model = data_structure.Model()
        coords = np.array([[0, 0, 0], [1, 0, 0]])
        mp_to = model.create_model_part('mp_to_1', *split(coords), np.arange(2))

        coords = np.vstack((coords, np.array([[1e-10, 0., 0.]])))
        mp_from = model.create_model_part('mp_from_1', *split(coords), np.arange(3))
        self.assertRaises(Warning, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()

        coords = np.vstack((coords, np.array([[1e-14, 0., 0.]])))
        mp_from = model.create_model_part('mp_from_2', *split(coords), np.arange(4))
        self.assertRaises(ValueError, mapper.initialize, *(mp_from, mp_to))
        mapper.finalize()


class Case1D:
    # 1D case: square-root grid + linear function
    def __init__(self, n_from, n_to):
        self.n_from = n_from
        self.n_to = n_to
        self.var = 'pressure'
        self.mp_name_from = 'wall_from'
        self.mp_name_to = 'wall_to'

        self.model = data_structure.Model()

        # Interface from
        self.z_from = np.linspace(0, 10, self.n_from) ** .5
        self.v_from = self.fun(self.z_from)
        self.model.create_model_part(self.mp_name_from, np.zeros(self.n_from),
                                     np.zeros(self.n_from), self.z_from, np.arange(self.n_from))
        parameters = [{'model_part': self.mp_name_from, 'variables': [self.var]}]
        self.interface_from = data_structure.Interface(parameters, self.model)
        self.interface_from.set_interface_data(self.v_from)

        # Interface to
        self.z_to = np.linspace(0, 10, self.n_to) ** .5
        self.model.create_model_part(self.mp_name_to, np.zeros(self.n_to),
                                     np.zeros(self.n_to), self.z_to, np.arange(self.n_to))
        parameters = [{'model_part': self.mp_name_to, 'variables': [self.var]}]
        self.interface_to = data_structure.Interface(parameters, self.model)

    def map(self, parameters):
        mapper = create_instance(parameters)
        mapper.initialize(self.model.get_model_part(self.mp_name_from),
                          self.model.get_model_part(self.mp_name_to))

        args_from = (self.interface_from, self.mp_name_from, self.var)
        args_to = (self.interface_to, self.mp_name_to, self.var)
        mapper(args_from, args_to)

        self.v_to_fun = self.fun(self.z_to)
        self.v_to = self.interface_to.get_variable_data(self.mp_name_to, self.var)

        self.v_error = np.abs(self.v_to.flatten() - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))
        plt.title(f'max error = {self.v_error.max():.2g}')
        plt.suptitle('1D case: square-root grid + linear function')

        ax[0].plot(self.z_from, self.v_from, label='from', marker='o')
        ax[0].plot(self.z_to, self.v_to, label='to', marker='o')
        ax[1].plot(self.z_to, self.v_error, label='error', marker='o')
        for a in ax:
            a.legend()
            a.set_xlabel('z')
            a.set_ylabel('f(z)')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, z):
        return z / 10


class Case2D:
    # 2D case: circle + linear function
    def __init__(self, n_from, n_to):
        self.n_from = n_from
        self.n_to = n_to
        self.var = 'pressure'
        self.mp_name_from = 'wall_from'
        self.mp_name_to = 'wall_to'

        self.model = data_structure.Model()

        # Interface from
        dtheta = 2 * np.pi / self.n_from
        self.theta_from = np.linspace(0, 2 * np.pi - dtheta, self.n_from)
        self.x_from, self.y_from = self.get_cartesian(self.theta_from)
        self.v_from = self.fun(self.x_from, self.y_from)
        self.model.create_model_part(self.mp_name_from, self.x_from, self.y_from,
                                     np.zeros(self.n_from), np.arange(self.n_from))
        parameters = [{'model_part': self.mp_name_from, 'variables': [self.var]}]
        self.interface_from = data_structure.Interface(parameters, self.model)
        self.interface_from.set_interface_data(self.v_from)

        # Interface to
        dtheta = 2 * np.pi / self.n_to
        self.theta_to = np.linspace(0, 2 * np.pi - dtheta, self.n_to)
        self.x_to, self.y_to = self.get_cartesian(self.theta_to)
        self.model.create_model_part(self.mp_name_to, self.x_to, self.y_to,
                                     np.zeros(self.n_to), np.arange(self.n_to))
        parameters = [{'model_part': self.mp_name_to, 'variables': [self.var]}]
        self.interface_to = data_structure.Interface(parameters, self.model)

    def map(self, parameters):
        mapper = create_instance(parameters)
        mapper.initialize(self.model.get_model_part(self.mp_name_from),
                          self.model.get_model_part(self.mp_name_to))

        args_from = (self.interface_from, self.mp_name_from, self.var)
        args_to = (self.interface_to, self.mp_name_to, self.var)
        mapper(args_from, args_to)

        self.v_to_fun = self.fun(self.x_to, self.y_to)
        self.v_to = self.interface_to.get_variable_data(self.mp_name_to, self.var)

        self.v_error = np.abs(self.v_to.flatten() - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        _, ax = plt.subplots(ncols=2, sharex=True, figsize=(15, 6))
        plt.suptitle('2D case: circle + linear function')

        ax[0].plot(self.theta_from * 180 / np.pi, self.v_from, label='from', marker='o')
        ax[0].plot(self.theta_to * 180 / np.pi, self.v_to, label='to', marker='o')

        ax[1].plot(self.theta_to * 180 / np.pi, self.v_error, label='error', marker='o')

        for a in ax:
            a.legend()
            a.set_xlabel(r'$\theta$')
            a.set_ylabel(r'f($\theta$)')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, x, y):
        return 2 * x + 3 * y

    def get_cartesian(self, theta):
        r = 2.
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return x, y


class Case3DSphere:
    # 3D case: sphere + sine function
    def __init__(self, n_theta_from, n_phi_from, n_theta_to, n_phi_to):
        self.n_theta_from = n_theta_from
        self.n_phi_from = n_phi_from
        self.n_from = n_theta_from * n_phi_from
        self.n_theta_to = n_theta_to
        self.n_phi_to = n_phi_to  # for bounding box: not too far from n_phi_from!
        self.n_to = n_theta_to * n_phi_to
        self.var = 'pressure'
        self.mp_name_from = 'wall_from'
        self.mp_name_to = 'wall_to'
        self.model = data_structure.Model()

        # Interface from
        shape = (self.n_theta_from, self.n_phi_from)
        dtheta = 2 * np.pi / self.n_theta_from
        dphi = np.pi / (self.n_phi_from - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_from).reshape(1, -1)

        self.x_from, self.y_from, self.z_from = self.get_cartesian(theta, phi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        self.model.create_model_part(self.mp_name_from, self.x_from.flatten(), self.y_from.flatten(),
                                     self.z_from.flatten(), np.arange(self.n_from))
        parameters = [{'model_part': self.mp_name_from, 'variables': [self.var]}]
        self.interface_from = data_structure.Interface(parameters, self.model)
        self.interface_from.set_interface_data(self.v_from.flatten())

        # Interface to
        shape = (self.n_theta_to, self.n_phi_to)
        dtheta = 2 * np.pi / self.n_theta_to
        dphi = np.pi / (self.n_phi_to - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_to).reshape(1, -1)

        self.x_to, self.y_to, self.z_to = self.get_cartesian(theta, phi)
        self.model.create_model_part(self.mp_name_to, self.x_to.flatten(), self.y_to.flatten(),
                                     self.z_to.flatten(), np.arange(self.n_to))
        parameters = [{'model_part': self.mp_name_to, 'variables': [self.var]}]
        self.interface_to = data_structure.Interface(parameters, self.model)

    def map(self, parameters):
        mapper = create_instance(parameters)
        mapper.initialize(self.model.get_model_part(self.mp_name_from),
                          self.model.get_model_part(self.mp_name_to))

        args_from = (self.interface_from, self.mp_name_from, self.var)
        args_to = (self.interface_to, self.mp_name_to, self.var)
        mapper(args_from, args_to)

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        self.v_to = self.interface_to.get_variable_data(self.mp_name_to, self.var)
        self.v_to = self.v_to.reshape(self.x_to.shape)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

    def check(self, tolerance):
        criterion = (self.v_error < tolerance)
        return criterion.all()

    def plot(self):
        v_min = min(self.v_from.min(), self.v_to.min())
        v_max = max(self.v_from.max(), self.v_to.max())
        c_from = cm.jet((self.v_from - v_min) / (v_max - v_min))
        c_to = cm.jet((self.v_to - v_min) / (v_max - v_min))
        c_error = cm.jet(self.v_error / self.v_error.max())

        fig = plt.figure(figsize=(18, 6))
        plt.suptitle(f'3D case: sphere + sine function | max error = {self.v_error.max():.2g}     ({v_min:.1f} < v < {v_max:.1g})')

        ax_from = fig.add_subplot(131, projection='3d')
        ax_from.set_title('from')
        ax_from.plot_surface(self.x_from, self.y_from, self.z_from, facecolors=c_from,
                             rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_to = fig.add_subplot(132, projection='3d')
        ax_to.set_title('to')
        ax_to.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_to,
                           rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_error = fig.add_subplot(133, projection='3d')
        ax_error.set_title('to (error)')
        ax_error.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_error,
                              rstride=1, cstride=1, antialiased=False, shade=False)

        for ax in [ax_from, ax_to, ax_error]:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        plt.tight_layout()
        plt.show()
        plt.close()

    def fun(self, x, y, z):
        return np.sin(x) * np.sin(y) * np.sin(z)

    def get_cartesian(self, theta, phi):
        r = np.pi
        x = r * np.cos(theta) * np.sin(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(phi)
        return x, y, z


class Case3DCylinder(Case3DSphere):
    # 3D case: cylinder + sine function
    def __init__(self, n_x_from, n_theta_from, n_x_to, n_theta_to, length):
        self.n_x_from = n_x_from
        self.n_theta_from = n_theta_from
        self.n_from = n_x_from * n_theta_from
        self.n_x_to = n_x_to
        self.n_theta_to = n_theta_to
        self.n_to = n_x_to * n_theta_to
        self.length = length
        self.var = 'pressure'
        self.mp_name_from = 'wall_from'
        self.mp_name_to = 'wall_to'
        self.model = data_structure.Model()

        self.model = data_structure.Model()

        # Interface from
        shape = (self.n_x_from, self.n_theta_from)
        dtheta = 2 * np.pi / self.n_theta_from
        theta_from = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(1, -1)

        self.x_from = np.ones(shape) * np.linspace(0, self.length, self.n_x_from).reshape(-1, 1)
        self.y_from, self.z_from = self.get_cartesian(theta_from)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        self.model.create_model_part(self.mp_name_from, self.x_from.flatten(), self.y_from.flatten(),
                                     self.z_from.flatten(), np.arange(self.n_from))
        parameters = [{'model_part': self.mp_name_from, 'variables': [self.var]}]
        self.interface_from = data_structure.Interface(parameters, self.model)
        self.interface_from.set_interface_data(self.v_from.flatten())

        # Interface to
        shape = (self.n_x_to, self.n_theta_to)
        dtheta = 2 * np.pi / self.n_theta_to
        theta_to = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(1, -1)

        self.x_to = np.ones(shape) * np.linspace(0, self.length, self.n_x_to).reshape(-1, 1)
        self.y_to, self.z_to = self.get_cartesian(theta_to)
        self.model.create_model_part(self.mp_name_to, self.x_to.flatten(), self.y_to.flatten(),
                                     self.z_to.flatten(), np.arange(self.n_to))
        parameters = [{'model_part': self.mp_name_to, 'variables': [self.var]}]
        self.interface_to = data_structure.Interface(parameters, self.model)

    def plot(self):
        v_min = min(self.v_from.min(), self.v_to.min())
        v_max = max(self.v_from.max(), self.v_to.max())
        c_from = cm.jet((self.v_from - v_min) / (v_max - v_min))
        c_to = cm.jet((self.v_to - v_min) / (v_max - v_min))
        c_error = cm.jet(self.v_error / self.v_error.max())

        fig = plt.figure(figsize=(18, 10))
        plt.suptitle(f'3D case: cylinder + sine function | max error = {self.v_error.max():.2g}     ({v_min:.1f} < v < {v_max:.1g})')

        # 2D plots
        ax_2dval = fig.add_subplot(221)
        ax_2dval.plot(self.x_from[:, 0], self.v_from[:, 0], label='from', marker='o')
        ax_2dval.plot(self.x_to[:, 0], self.v_to[:, 0], label='to', marker='o')

        ax_2derr = fig.add_subplot(222)
        ax_2derr.plot(self.x_to[:, 0], self.v_error[:, 0], label='error', marker='o')

        for ax in [ax_2dval, ax_2derr]:
            ax.legend()
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'f($x$)')

        # 3D plots
        ax_from = fig.add_subplot(234, projection='3d')
        ax_from.set_title('from')
        ax_from.plot_surface(self.x_from, self.y_from, self.z_from, facecolors=c_from,
                             rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_to = fig.add_subplot(235, projection='3d')
        ax_to.set_title('to')
        ax_to.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_to,
                           rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

        ax_error = fig.add_subplot(236, projection='3d')
        ax_error.set_title('to (error)')
        ax_error.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_error,
                              rstride=1, cstride=1, antialiased=False, shade=False)

        for ax in [ax_from, ax_to, ax_error]:
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        plt.show()
        plt.close()

    def fun(self, x, y, z):
        theta = np.arctan2(z, y)
        return np.sin(2 * np.pi * x / self.length) * np.cos(theta)

    def get_cartesian(self, theta):
        r = .5
        y = r * np.cos(theta)
        z = r * np.sin(theta)
        return y, z


class Case3DSinc:
    # 3D case: sinc + linear vector function
    def __init__(self, n_x_from, n_y_from, n_x_to, n_y_to):
        self.n_x_from = n_x_from
        self.n_y_from = n_y_from
        self.n_from = n_x_from * n_y_from
        self.n_x_to = n_x_to
        self.n_y_to = n_y_to
        self.n_to = n_x_to * n_y_to
        self.var = 'displacement'
        self.mp_name_from = 'wall_from'
        self.mp_name_to = 'wall_to'
        self.model = data_structure.Model()

        model = data_structure.Model()

        # ModelPart from
        shape = (self.n_x_from, self.n_y_from)
        self.x_from = np.ones(shape) * np.linspace(-1, 1, self.n_x_from).reshape(-1, 1)
        self.y_from = np.ones(shape) * np.linspace(-1, 1, self.n_y_from).reshape(1, -1)
        self.z_from = np.sinc(np.sqrt((10 * self.x_from) ** 2 + (self.y_from * 10) ** 2) / np.pi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        self.model.create_model_part(self.mp_name_from, self.x_from.flatten(), self.y_from.flatten(),
                                     self.z_from.flatten(), np.arange(self.n_from))
        parameters = [{'model_part': self.mp_name_from, 'variables': [self.var]}]
        self.interface_from = data_structure.Interface(parameters, self.model)
        tmp = np.hstack((self.v_from[0].reshape(-1, 1),
                         self.v_from[1].reshape(-1, 1),
                         self.v_from[2].reshape(-1, 1)))
        self.interface_from.set_variable_data(self.mp_name_from, self.var, tmp)

        # ModelPart to
        shape = (self.n_x_to, self.n_y_to)
        self.x_to = np.ones(shape) * np.linspace(-1, 1, self.n_x_to).reshape(-1, 1)
        self.y_to = np.ones(shape) * np.linspace(-1, 1, self.n_y_to).reshape(1, -1)
        self.z_to = np.sinc(np.sqrt((10 * self.x_to) ** 2 + (self.y_to * 10) ** 2) / np.pi)
        self.v_to = self.fun(self.x_to, self.y_to, self.z_to)
        self.model.create_model_part(self.mp_name_to, self.x_to.flatten(), self.y_to.flatten(),
                                     self.z_to.flatten(), np.arange(self.n_to))
        parameters = [{'model_part': self.mp_name_to, 'variables': [self.var]}]
        self.interface_to = data_structure.Interface(parameters, self.model)

    def map(self, parameters):
        mapper = create_instance(parameters)
        mapper.initialize(self.model.get_model_part(self.mp_name_from),
                          self.model.get_model_part(self.mp_name_to))

        args_from = (self.interface_from, self.mp_name_from, self.var)
        args_to = (self.interface_to, self.mp_name_to, self.var)
        mapper(args_from, args_to)

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        tmp = self.interface_to.get_variable_data(self.mp_name_to, self.var)
        shape_to = self.v_to_fun[0].shape
        self.v_to = [tmp[:, 0].reshape(shape_to),
                     tmp[:, 1].reshape(shape_to),
                     tmp[:, 2].reshape(shape_to)]

        self.v_error = []
        for j in range(3):
            self.v_error.append(np.abs(self.v_to[j] - self.v_to_fun[j]))

    def check(self, tolerance):
        out = []
        for j in range(3):
            criterion = (self.v_error[j] < tolerance)
            out.append(criterion.all())
        return out

    def plot(self):
        fig = plt.figure(figsize=(18, 10))
        plt.suptitle('3D case: sinc + linear vector function')

        for j in range(3):
            v_min = min(self.v_from[j].min(), self.v_to[j].min())
            v_max = max(self.v_from[j].max(), self.v_to[j].max())
            c_from = cm.jet((self.v_from[j] - v_min) / (v_max - v_min))
            c_to = cm.jet((self.v_to[j] - v_min) / (v_max - v_min))
            c_error = cm.jet(self.v_error[j] / self.v_error[j].max())

            ax_from = fig.add_subplot(3, 3, 1 + j * 3, projection='3d')
            ax_from.set_title('from')
            ax_from.plot_surface(self.x_from, self.y_from, self.z_from, facecolors=c_from,
                                 rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

            ax_to = fig.add_subplot(3, 3, 2 + j * 3, projection='3d')
            ax_to.set_title('to')
            ax_to.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_to,
                               rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)

            ax_error = fig.add_subplot(3, 3, 3 + j * 3, projection='3d')
            ax_error.set_title(f'max error = {self.v_error[j].max():.2g} ({v_min:.1f} < v < {v_max:.1g})')
            ax_error.plot_surface(self.x_to, self.y_to, self.z_to, facecolors=c_error,
                                  rstride=1, cstride=1, antialiased=False, shade=False)

            for ax in [ax_from, ax_to, ax_error]:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')

        plt.tight_layout()
        plt.get_current_fig_manager().window.showMaximized()
        plt.show()
        plt.close()

    def fun(self, x, y, z):
        return [y + z, z + x, x + y]


if __name__ == '__main__':
    unittest.main()
