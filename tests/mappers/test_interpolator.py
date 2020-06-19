from coconut import data_structure
from coconut.data_structure import KratosUnittest
from coconut.coupling_components.tools import CreateInstance

import os
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

variables = vars(data_structure)


class TestMapperInterpolator(KratosUnittest.TestCase):
    def test_mapper_interpolator(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_interpolator.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = data_structure.Parameters(parameter_file.read())
        par_mapper_0 = parameters['mapper']

        # test if directions are set correctly in __init__
        par_mapper = deepcopy(par_mapper_0)

        par_mapper['settings'].SetArray('directions', ['Z'])
        mapper = CreateInstance(par_mapper)
        self.assertListEqual(mapper.directions, ['Z0'])

        par_mapper['settings'].SetArray('directions', ['Z', 'Y', 'X'])
        mapper = CreateInstance(par_mapper)
        self.assertListEqual(mapper.directions, ['Z0', 'Y0', 'X0'])

        par_mapper['settings'].SetArray('directions', ['z'])
        self.assertRaises(ValueError, CreateInstance, par_mapper)

        par_mapper['settings'].SetArray('directions', ['Z0'])
        self.assertRaises(ValueError, CreateInstance, par_mapper)

        par_mapper['settings'].SetArray('directions', ['Z', 'Y', 'X', 'Z'])
        self.assertRaises(ValueError, CreateInstance, par_mapper)

        par_mapper['settings'].SetString('directions', 'Z')
        self.assertRaises(TypeError, CreateInstance, par_mapper)

        # test check_bounding_box method
        if True:
            par_mapper = deepcopy(par_mapper_0)

            # check 1D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z'])
            mapper = CreateInstance(par_mapper)

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 0., 0., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 0., 0., 1.)

            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(2, 0., 0., 1.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(3, 0., 0., -.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(11, 0., 0., 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(12, 0., 0., 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(13, 0., 0., -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(14, 0., 0., 2.)
            mp_to.CreateNewNode(15, 0., 0., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            # check 2D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z', 'X'])
            mapper = CreateInstance(par_mapper)

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 1., 0., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 1., 0., 1.)

            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(2, 1.01, 0., 1.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(3, -.01, 0., -.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(11, 1.1, 0., 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(12, 1.25, 0., 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(13, -.25, 0., -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(14, 2., 0., 2.)
            mp_to.CreateNewNode(15, -1., 0., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            # check 3D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z', 'X', 'Y'])
            mapper = CreateInstance(par_mapper)

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 1., 1., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 1., 1., 1.)

            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(2, 1.01, 1.01, 1.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(3, -.01, -.01, -.01)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            mp_to.CreateNewNode(11, 1.1, 1.1, 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(12, 1.25, 1.25, 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(13, -.25, -.25, -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            mp_to.CreateNewNode(14, 2., 2., 2.)
            mp_to.CreateNewNode(15, -1., -1., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            # check if method works for lines aligned with coordinate axes in 2D
            par_mapper['settings'].SetArray('directions', ['X', 'Y'])
            mapper = CreateInstance(par_mapper)

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1., 0.)
            mp_to.CreateNewNode(1, 1.01, 1.01, 0.)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.05, 0.)
            mp_to.CreateNewNode(1, 1., 1.05, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.25, 0.)
            mp_to.CreateNewNode(1, 1., 1.25, 0.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            # mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., .85, 0.)
            mp_to.CreateNewNode(1, 1., 1.15, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            # check if method works for planes aligned with coordinate axes in 3D
            par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
            mapper = CreateInstance(par_mapper)

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1., 0.)
            mp_to.CreateNewNode(1, 1.01, 1.01, 0.)
            mapper.Initialize(mp_from, mp_to)
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.05, 0.)
            mp_to.CreateNewNode(1, 1., 1.05, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.25, 0.)
            mp_to.CreateNewNode(1, 1., 1.25, 0.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

            model = data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., .85, 0.)
            mp_to.CreateNewNode(1, 1., 1.15, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
            mapper.Finalize()

        # test check_duplicate_points method
        par_mapper = deepcopy(par_mapper_0)
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        mapper = CreateInstance(par_mapper)

        model = data_structure.Model()
        mp_from = model.CreateModelPart('mp_from')
        mp_to = model.CreateModelPart('mp_to')
        mp_from.CreateNewNode(0, 0., 0., 0.)
        mp_from.CreateNewNode(1, 1., 0., 0.)
        mp_to.CreateNewNode(0, 0., 0., 0.)
        mp_to.CreateNewNode(1, 1., 0., 0.)

        mp_from.CreateNewNode(2, 1e-10, 0., 0.)
        self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))
        mapper.Finalize()

        mp_from.CreateNewNode(3, 1e-14, 0., 0.)
        self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
        mapper.Finalize()

        # to do: check tree? check __call__ method?


class Case1D:
    # 1D case: square-root grid + linear function
    def __init__(self, n_from, n_to):
        self.n_from = n_from
        self.n_to = n_to

        model = data_structure.Model()

        # ModelPart from
        self.var_from = variables["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)
        self.z_from = np.linspace(0, 10, self.n_from) ** .5
        self.v_from = self.fun(self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, 0., 0., self.z_from[i])
            node.SetSolutionStepValue(self.var_from, 0, self.v_from[i])

        # ModelPart to
        self.var_to = variables["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)
        self.z_to = np.linspace(0, 10, self.n_to) ** .5
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, 0., 0., self.z_to[i])

    def map(self, parameters):
        mapper = CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.z_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

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

        model = data_structure.Model()

        # ModelPart from
        self.var_from = variables["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        dtheta = 2 * np.pi / self.n_from
        self.theta_from = np.linspace(0, 2 * np.pi - dtheta, self.n_from)

        self.x_from, self.y_from = self.get_cartesian(self.theta_from)
        self.v_from = self.fun(self.x_from, self.y_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from[i], self.y_from[i], 0.)
            node.SetSolutionStepValue(self.var_from, 0, self.v_from[i])

        # ModelPart to
        self.var_to = variables["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        dtheta = 2 * np.pi / self.n_to
        self.theta_to = np.linspace(0, 2 * np.pi - dtheta, self.n_to)

        self.x_to, self.y_to = self.get_cartesian(self.theta_to)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to[i], self.y_to[i], 0.)

    def map(self, parameters):
        mapper = CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
        self.v_error = np.abs(self.v_to - self.v_to_fun)

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

        model = data_structure.Model()

        # ModelPart from
        self.var_from = variables["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_theta_from, self.n_phi_from)
        dtheta = 2 * np.pi / self.n_theta_from
        dphi = np.pi / (self.n_phi_from - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_from).reshape(1, -1)

        self.x_from, self.y_from, self.z_from = self.get_cartesian(theta, phi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from.flatten()[i],
                                self.y_from.flatten()[i], self.z_from.flatten()[i])
            node.SetSolutionStepValue(self.var_from, 0, self.v_from.flatten()[i])

        # ModelPart to
        self.var_to = variables["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        shape = (self.n_theta_to, self.n_phi_to)
        dtheta = 2 * np.pi / self.n_theta_to
        dphi = np.pi / (self.n_phi_to - 1)
        theta = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(-1, 1)
        phi = np.ones(shape) * np.linspace(dphi, np.pi - dphi, self.n_phi_to).reshape(1, -1)

        self.x_to, self.y_to, self.z_to = self.get_cartesian(theta, phi)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                            self.y_to.flatten()[i], self.z_to.flatten()[i])

    def map(self, parameters):
        mapper = CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        self.v_to = np.zeros(self.n_to)
        for i, node in enumerate(self.model_part_to.Nodes):
            self.v_to[i] = node.GetSolutionStepValue(self.var_to)
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

        model = data_structure.Model()

        # ModelPart from
        self.var_from = variables["TEMPERATURE"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_x_from, self.n_theta_from)
        dtheta = 2 * np.pi / self.n_theta_from
        theta_from = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_from).reshape(1, -1)

        self.x_from = np.ones(shape) * np.linspace(0, self.length, self.n_x_from).reshape(-1, 1)
        self.y_from, self.z_from = self.get_cartesian(theta_from)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from.flatten()[i],
                                self.y_from.flatten()[i], self.z_from.flatten()[i])
            node.SetSolutionStepValue(self.var_from, 0, self.v_from.flatten()[i])
        # for i in range(self.n_from):
        #     node = self.model_part_from.CreateNewNode(i, self.x_from[i], self.y_from[i], self.z_from[i])
        #     node.SetSolutionStepValue(self.var_from, 0, self.v_from[i])

        # ModelPart to
        self.var_to = variables["PRESSURE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        shape = (self.n_x_to, self.n_theta_to)
        dtheta = 2 * np.pi / self.n_theta_to
        theta_to = np.ones(shape) * np.linspace(0, 2 * np.pi - dtheta, self.n_theta_to).reshape(1, -1)

        self.x_to = np.ones(shape) * np.linspace(0, self.length, self.n_x_to).reshape(-1, 1)
        self.y_to, self.z_to = self.get_cartesian(theta_to)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                            self.y_to.flatten()[i], self.z_to.flatten()[i])
        # for i in range(self.n_to):
        #     self.model_part_to.CreateNewNode(i, self.x_to[i], self.y_to[i], self.z_to[i])

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

        plt.tight_layout()
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

        model = data_structure.Model()

        # ModelPart from
        self.var_from = variables["VELOCITY"]
        self.model_part_from = model.CreateModelPart('wall_from')
        self.model_part_from.AddNodalSolutionStepVariable(self.var_from)

        shape = (self.n_x_from, self.n_y_from)
        self.x_from = np.ones(shape) * np.linspace(-1, 1, self.n_x_from).reshape(-1, 1)
        self.y_from = np.ones(shape) * np.linspace(-1, 1, self.n_y_from).reshape(1, -1)
        self.z_from = np.sinc(np.sqrt((10 * self.x_from) ** 2 + (self.y_from * 10) ** 2) / np.pi)
        self.v_from = self.fun(self.x_from, self.y_from, self.z_from)
        for i in range(self.n_from):
            node = self.model_part_from.CreateNewNode(i, self.x_from.flatten()[i],
                                    self.y_from.flatten()[i], self.z_from.flatten()[i])
            hist = []
            for j in range(3):
                hist.append(self.v_from[j].flatten()[i])
            node.SetSolutionStepValue(self.var_from, 0, hist)

        # ModelPart to
        self.var_to = variables["FORCE"]
        self.model_part_to = model.CreateModelPart('wall_to')
        self.model_part_to.AddNodalSolutionStepVariable(self.var_to)

        shape = (self.n_x_to, self.n_y_to)
        self.x_to = np.ones(shape) * np.linspace(-1, 1, self.n_x_to).reshape(-1, 1)
        self.y_to = np.ones(shape) * np.linspace(-1, 1, self.n_y_to).reshape(1, -1)
        self.z_to = np.sinc(np.sqrt((10 * self.x_to) ** 2 + (self.y_to * 10) ** 2) / np.pi)
        self.v_to = self.fun(self.x_to, self.y_to, self.z_to)
        for i in range(self.n_to):
            self.model_part_to.CreateNewNode(i, self.x_to.flatten()[i],
                                 self.y_to.flatten()[i], self.z_to.flatten()[i])

    def map(self, parameters):
        mapper = CreateInstance(parameters)
        mapper.Initialize(self.model_part_from, self.model_part_to)
        mapper((self.model_part_from, self.var_from),
               (self.model_part_to, self.var_to))

        self.v_to_fun = self.fun(self.x_to, self.y_to, self.z_to)
        self.v_to = [np.zeros(self.n_to), np.zeros(self.n_to), np.zeros(self.n_to)]
        for i, node in enumerate(self.model_part_to.Nodes):
            hist = node.GetSolutionStepValue(self.var_to)
            for j in range(3):
                self.v_to[j][i] = hist[j]

        self.v_error = []
        for j in range(3):
            self.v_to[j] = self.v_to[j].reshape(self.x_to.shape)
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
    KratosUnittest.main()