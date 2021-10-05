from coconut import data_structure
from coconut.tools import create_instance

import numpy as np
import unittest
import matplotlib.pyplot as plt
from matplotlib import cm


class TestMapperAxisymmetric2DTo3D(unittest.TestCase):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.axisymmetric_2d_to_3d',
                           'settings':
                               {'direction_axial': 'x',
                                'direction_radial': 'y',
                                'angle': 100,
                                'n_tangential': 8}
                           }
        self.forward = True

    def test_instantiation(self):
        create_instance(self.parameters)

        self.parameters['settings']['direction_axial'] = 'X'
        self.assertRaises(ValueError, create_instance, self.parameters)

        self.parameters['settings']['direction_axial'] = 'x'
        self.parameters['settings']['direction_radial'] = 2
        self.assertRaises(ValueError, create_instance, self.parameters)

        self.parameters['settings']['direction_radial'] = 2
        self.parameters['settings']['n_tangential'] = 3
        self.assertRaises(ValueError, create_instance, self.parameters)

    def test_initialize(self):
        mp_name_in = 'wall_in'
        mp_name_out = 'wall_out'

        # create model_part_in
        n_in = 10
        x_in = np.linspace(0, 2 * np.pi, n_in)
        y_in = 1. + 0.2 * np.sin(x_in)
        z_in = np.zeros(n_in)
        model = data_structure.Model()
        model.create_model_part(mp_name_in, x_in, y_in, z_in, np.arange(n_in))

        # initialize mapper to get model_part_out
        mapper = create_instance(self.parameters)
        with self.assertRaises(NotImplementedError):
            mapper.initialize(model, mp_name_in, mp_name_out, forward=not self.forward)
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out, forward=self.forward)

        # get mapped geometry from 3D model_part_out
        mp_out = model.get_model_part(mp_name_out)
        x_out = mp_out.x0
        y_out = mp_out.y0
        z_out = mp_out.z0

        # check mapped geometry
        n_t = self.parameters['settings']['n_tangential']
        total_angle = self.parameters['settings'].get('angle', 360)
        for x_i in x_in:
            y = y_out[x_out == x_i]
            z = z_out[x_out == x_i]

            # check radius
            radius_in = np.sqrt(y_in[x_in == x_i] ** 2 + z_in[x_in == x_i] ** 2)
            radii_in = np.full((n_t,), radius_in)
            radii_out = np.sqrt(y ** 2 + z ** 2)
            np.testing.assert_allclose(radii_in, radii_out, rtol=1e-14)

            # check angle between consecutive points
            angle_in = total_angle / n_t if total_angle == 360 else total_angle / (n_t - 1)
            angles_in = np.full((n_t - 1,), angle_in)
            theta = np.arccos(y / radius_in) * np.sign(z) * 180 / np.pi
            np.sort(theta)
            angles_out = np.diff(theta)
            np.testing.assert_allclose(angles_in, angles_out, rtol=1e-14)

            # check total (wedge) angle of new model part
            total_angle_in = total_angle - angle_in if total_angle == 360 else total_angle
            total_angle_out = theta[-1] - theta[0]
            self.assertAlmostEqual(total_angle_out, total_angle_in)

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

        n_from = 10
        tmp = np.linspace(0, 5, n_from)
        x_from, y_from, z_from = tmp, 1. + 0.2 * np.sin(2 * np.pi / 5 * tmp), np.zeros_like(tmp)
        v_s_from = fun_s(x_from).reshape(-1, 1)
        v_v_from = fun_v(x_from, y_from, z_from)

        model = data_structure.Model()
        model.create_model_part(mp_name_from, x_from, y_from, y_from, np.arange(n_from))
        parameters_from = [{'model_part': mp_name_from, 'variables': [var_s, var_v]}]
        interface_from = data_structure.Interface(parameters_from, model)
        interface_from.set_variable_data(mp_name_from, var_s, v_s_from)
        interface_from.set_variable_data(mp_name_from, var_v, v_v_from)

        # initialize mapper
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_from, mp_name_to, forward=True)
        parameters_to = [{'model_part': mp_name_to, 'variables': [var_s, var_v]}]
        interface_to = data_structure.Interface(parameters_to, model)
        mp_to = interface_to.get_model_part(mp_name_to)
        x_to, y_to, z_to = mp_to.x0, mp_to.y0, mp_to.z0

        # check mapped values for 1D variable
        mapper((interface_from, mp_name_from, var_s),
               (interface_to, mp_name_to, var_s))
        v_s_to_ref = fun_s(x_to).reshape(-1, 1)
        v_s_to = interface_to.get_variable_data(mp_name_to, var_s)
        np.testing.assert_allclose(v_s_to, v_s_to_ref, rtol=1e-14)

        # check mapped values for 3D variable
        mapper((interface_from, mp_name_from, var_v),
               (interface_to, mp_name_to, var_v))
        v_v_to_ref = fun_v(x_to, y_to, z_to)
        v_v_to = interface_to.get_variable_data(mp_name_to, var_v)
        np.testing.assert_allclose(v_v_to, v_v_to_ref, rtol=1e-14)

        # extra: visualization
        if self.gui:
            v_s_from, v_s_to = v_s_from.flatten(), v_s_to.flatten()
            c_from = cm.jet((v_s_from - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))
            c_to = cm.jet((v_s_to - v_s_from.min()) / (v_s_from.max() - v_s_from.min()))

            fig = plt.figure()

            ax_s = fig.add_subplot(121, projection='3d')
            ax_s.set_title('check geometry and scalar mapping')
            ax_s.scatter(x_from, y_from, z_from, s=50, c=c_from, depthshade=True, marker='s')
            ax_s.scatter(x_to, y_to, z_to, s=20, c=c_to, depthshade=True)

            ax_v = fig.add_subplot(122, projection='3d')
            ax_v.set_title('check vector mapping')
            ax_v.quiver(x_from, y_from, z_from, v_v_from[:, 0], v_v_from[:, 1], v_v_from[:, 2],
                        pivot='tail', arrow_length_ratio=0.1, normalize=False, length=0.1, colors='r', linewidth=3)
            ax_v.quiver(x_to, y_to, z_to, v_v_to[:, 0], v_v_to[:, 1], v_v_to[:, 2],
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


if __name__ == '__main__':
    unittest.main()
