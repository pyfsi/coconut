import unittest

import matplotlib.pyplot as plt
import numpy as np
from coconut import data_structure
from coconut.tools import create_instance
from matplotlib import cm


class TestMapperAxisymmetric3DTo2D(unittest.TestCase):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.axisymmetric_3d_to_2d',
                           'settings':
                               {'direction_axial': 'x',
                                'direction_radial': 'y',
                                'angle': 180,
                                'n_tangential': 9}
                           }

    def test_initialize(self):
        mp_name_in = 'wall_in'
        mp_name_out = 'wall_out'
        angle = self.parameters['settings'].get('angle', 360)

        # create model_part_in
        n_in = 10
        x_in = np.linspace(0, 2 * np.pi, n_in)
        y_in = 1. + 0.2 * np.sin(x_in)
        z_in = np.zeros(10)
        model = data_structure.Model()
        model.create_model_part(mp_name_in, x_in, y_in, z_in, np.arange(n_in))

        # create reference geometry for 3D model_part_out

        n_t = self.parameters['settings']['n_tangential']
        n_out_ref = n_in * n_t
        x_out_ref = np.zeros(n_out_ref)
        y_out_ref = np.zeros(n_out_ref)
        z_out_ref = np.zeros(n_out_ref)

        for i_t in range(n_t):
            for i_from in range(n_in):
                start = i_t * n_in
                end = (i_t + 1) * n_in
                if angle == 360:
                    theta = -np.radians(angle / 2) + i_t * np.radians(angle) / (n_t)
                else:
                    theta = -np.radians(angle / 2) + i_t * np.radians(angle) / (n_t - 1)
                x_out_ref[start:end] = x_in
                y_out_ref[start:end] = np.cos(theta) * y_in
                z_out_ref[start:end] = np.sin(theta) * y_in

        # test angle parameter
        radius = np.zeros(n_in)
        corner = np.zeros(n_in)

        for i_t in range(n_in):
            radius[i_t] = np.sqrt(y_out_ref[i_t] ** 2 + z_out_ref[i_t] ** 2)
            cosalpha = abs(y_out_ref[i_t]) / radius[i_t]
            corner[i_t] = 2 * np.arccos(cosalpha) * 180 / np.pi
            if angle > 180:
                corner[i_t] = 360 - corner[i_t]

        self.assertAlmostEqual(angle, corner[i_t], 10)

        # test angle between each point via cosine rule
        if angle != 360:
            ref = angle / (n_t - 1)
            a = np.zeros((n_t - 1) * n_in)
            b = np.zeros((n_t - 1) * n_in)
            c = np.zeros((n_t - 1) * n_in)

            k = 0

            for i_t in range(n_t - 1):
                for i in range(n_in):
                    a[k] = np.sqrt(y_out_ref[k] ** 2 + z_out_ref[k] ** 2)
                    b[k] = np.sqrt(y_out_ref[k + n_in] ** 2 + z_out_ref[k + n_in] ** 2)
                    c[k] = np.sqrt(
                        (y_out_ref[k] - y_out_ref[k + n_in]) ** 2 + (z_out_ref[k] - z_out_ref[k + n_in]) ** 2)
                    k += 1

            cos_gamma = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
            check_out = np.arccos(cos_gamma) * 180 / np.pi

            for i_t in range((n_t - 1) * n_in):
                self.assertAlmostEqual(ref, check_out[i_t], 10)

        else:
            ref = angle / n_t
            a = np.zeros(n_t * n_in)
            b = np.zeros(n_t * n_in)
            c = np.zeros(n_t * n_in)

            k = 0

            end_points = np.zeros(n_in)
            for i in range(n_in):
                end_points[i] = int(n_out_ref - n_in + i * 1)

            for i_t in range(n_t):
                for i in range(n_in):
                    a[k] = np.sqrt(y_out_ref[k] ** 2 + z_out_ref[k] ** 2)
                    if k in end_points:
                        b[k] = np.sqrt(
                            y_out_ref[k - ((n_t - 1) * n_in)] ** 2 + z_out_ref[k - ((n_t - 1) * n_in)] ** 2)
                        c[k] = np.sqrt((y_out_ref[k] - y_out_ref[k - ((n_t - 1) * n_in)]) ** 2 + (
                                z_out_ref[k] - z_out_ref[k - ((n_t - 1) * n_in)]) ** 2)
                    else:
                        b[k] = np.sqrt(y_out_ref[k + n_in] ** 2 + z_out_ref[k + n_in] ** 2)
                        c[k] = np.sqrt(
                            (y_out_ref[k] - y_out_ref[k + n_in]) ** 2 + (
                                    z_out_ref[k] - z_out_ref[k + n_in]) ** 2)
                    k += 1

            cos_gamma_360 = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
            check_out = np.arccos(cos_gamma_360) * 180 / np.pi

            for i_t in range(n_t * n_in):
                self.assertAlmostEqual(ref, check_out[i_t], 10)

        # initialize mapper to get model_part_out
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out, forward=False)

        # get mapped geometry from 3D model_part_out
        mp_out = model.get_model_part(mp_name_out)
        n_out = mp_out.size
        x_out = mp_out.x0
        y_out = mp_out.y0
        z_out = mp_out.z0

        # compare mapped and reference geometries
        self.assertEqual(n_out, n_out_ref)
        np.testing.assert_array_equal(x_out, x_out_ref)
        np.testing.assert_array_equal(y_out, y_out_ref)
        np.testing.assert_array_equal(z_out, z_out_ref)

    def test_call(self):
        def fun_s(x):
            return 1. + 2.5 * x

        def fun_v(x, y, z):
            theta = np.arctan2(z, y)
            v_x = 1. + 2.5 * x
            v_y = v_x * 0.5 * np.cos(theta)
            v_z = v_x * 0.5 * np.sin(theta)
            return np.column_stack((v_x, v_y, v_z))

        mp_name_from = 'wall_from'
        mp_name_to = 'wall_to'
        var_s = 'pressure'
        var_v = 'displacement'

        n_to = 10
        tmp = np.linspace(0, 5, n_to)
        x_to, y_to, z_to = tmp, 1. + 0.2 * np.sin(2 * np.pi / 5 * tmp), np.zeros_like(tmp)

        model = data_structure.Model()
        model.create_model_part(mp_name_to, x_to, y_to, z_to, np.arange(n_to))
        parameters_to = [{'model_part': mp_name_to, 'variables': [var_s, var_v]}]
        interface_to = data_structure.Interface(parameters_to, model)

        # initialize mapper
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_to, mp_name_from, forward=False)
        parameters_from = [{'model_part': mp_name_from, 'variables': [var_s, var_v]}]
        interface_from = data_structure.Interface(parameters_from, model)
        mp_from = interface_from.get_model_part(mp_name_from)
        x_from, y_from, z_from = mp_from.x0, mp_from.y0, mp_from.z0

        # set data
        v_s_from = fun_s(x_from).reshape(-1, 1)
        v_v_from = fun_v(x_from, y_from, z_from)
        interface_from.set_variable_data(mp_name_from, var_s, v_s_from)
        interface_from.set_variable_data(mp_name_from, var_v, v_v_from)

        # check mapped values for 1D variable
        mapper((interface_from, mp_name_from, var_s),
               (interface_to, mp_name_to, var_s))
        v_s_to_ref = fun_s(x_to).reshape(-1, 1)
        v_s_to = interface_to.get_variable_data(mp_name_to, var_s)
        np.testing.assert_allclose(v_s_to, v_s_to_ref, rtol=1e-15)

        # check mapped values for 3D variable
        mapper((interface_from, mp_name_from, var_v),
               (interface_to, mp_name_to, var_v))
        v_v_to_ref = fun_v(x_to, y_to, z_to)
        v_v_to = interface_to.get_variable_data(mp_name_to, var_v)
        np.testing.assert_allclose(v_v_to, v_v_to_ref, rtol=1e-15)

        # extra: visualization
        if self.gui:
            v_s_from, v_s_to = v_s_from.flatten(), v_s_to.flatten()
            c_from = cm.jet((v_s_from - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))
            c_to = cm.jet((v_s_to - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))

            fig = plt.figure()

            ax_s = fig.add_subplot(121, projection='3d')
            ax_s.set_title('check geometry and scalar mapping')
            ax_s.scatter(x_to, y_to, z_to, s=50, c=c_to, depthshade=True, marker='s')
            ax_s.scatter(x_from, y_from, z_from, s=20, c=c_from, depthshade=True)

            ax_v = fig.add_subplot(122, projection='3d')
            ax_v.set_title('check vector mapping')
            ax_v.quiver(x_to, y_to, z_to, v_v_to[:, 0], v_v_to[:, 1], v_v_to[:, 2],
                        pivot='tail', arrow_length_ratio=0.1, normalize=False, length=0.1, colors='r', linewidth=3)
            ax_v.quiver(x_from, y_from, z_from, v_v_from[:, 0], v_v_from[:, 1], v_v_from[:, 2],
                        pivot='tail', arrow_length_ratio=0.1, normalize=False, length=0.1)

            for ax in [ax_s, ax_v]:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')

            plt.get_current_fig_manager().window.showMaximized()
            plt.show()
            plt.close()

    def test_initialize_360(self):
        self.parameters['settings'].pop('angle')
        self.test_initialize()
        self.test_call()
        self.gui


if __name__ == '__main__':
    unittest.main()
