from coconut import data_structure
from coconut.tools import create_instance

import unittest
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import cm
from coconut.coupling_components.mappers.load_update2 import Mapper_load_update2


class TestMapper_load_update(unittest.TestCase):
    gui = True

    def setUp(self):
        self.parameters = {'type': 'mappers.load_update2',
                           'settings':
                               {"coords_min": 0.0,
                                "coords_max": 2.0}
                           }

    # def test_instantiation(self):
    #
    #     self.parameters['settings']["coords_min"] = 't'
    #     self.assertRaises(ValueError, Mapper_load_update2 , self.parameters)
    #
    #     self.parameters['settings']["coords_min"] = 0
    #     self.parameters['settings']["coords_max"] = 'y'
    #     self.assertRaises(ValueError, Mapper_load_update2, self.parameters)

    def test_initialize(self):
        x_min = -0.5
        x_max = 3
        mp_name_in = 'wall_in'
        mp_name_out = 'wall_out'


        # create extended geometry for model_part_in
        n_in = 10
        n_to = n_in
        x = np.linspace(x_min,x_max, n_in)
        r = 1 + 0.07 * np.sin(x * 600)

        x_in = np.zeros(n_to)
        y_in = np.zeros(n_to)
        z_in = np.zeros(n_to)

        i = 0
        for k in range(n_in):
            x_in[i] = x[k]
            y_in[i] = r[k]
            z_in[i] = 0
            i += 1

        model = data_structure.Model()
        self.mp_input_to= model.create_model_part(mp_name_in, x_in, y_in, z_in, np.arange(n_to))

        # initialize mapper to get model_part_out
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_in, mp_name_out, forward=False)

        # create new mp of size mp_input_to and length mp_input_from
        n_out_ref = self.mp_input_to.size
        x_out_ref = x_in
        y_out_ref = y_in
        z_out_ref = z_in

        # get mapped geometry from 2D model_part_out
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
            theta = np.arctan(y)
            v_x = 1. + 2.5 * x
            v_y = v_x * 0.5 * np.cos(theta)
            v_z = 0 * z
            return np.column_stack((v_x, v_y, v_z))

        mp_name_from = 'wall_from'
        mp_name_to = 'wall_to'
        var_s = 'pressure'
        var_v = 'displacement'

        n_in = 50
        n_to = n_in
        tmp = np.linspace(-0.5, 3, n_in)
        r_tmp = (1. + 0.2 * np.sin(2 * np.pi / 5 * tmp))

        # create model_part interface_input_to
        x_to = np.zeros(n_to)
        y_to = np.zeros(n_to)
        z_to = np.zeros(n_to)

        i = 0
        for k in range(n_to):
            x_to[i] = tmp[k]
            y_to[i] = r_tmp[k]
            z_to[i] = 0
            i += 1

        model = data_structure.Model()
        model.create_model_part(mp_name_to, x_to, y_to, z_to, np.arange(n_to))
        self.mp_to = model.get_model_part(mp_name_to)

        #initialize mapper to get mp_from
        mapper = create_instance(self.parameters)
        mapper.initialize(model, mp_name_to, mp_name_from, forward=False)

        #input_from interface simulated

        self.v_min = self.parameters['settings']['coords_min']
        self.v_max = self.parameters['settings']['coords_max']

        # setting of variables mp_from
        parameters_from = [{'model_part': mp_name_from, 'variables': [var_s, var_v]}]
        interface_from = data_structure.Interface(parameters_from, model)
        mp_from = interface_from.get_model_part(mp_name_from)
        x_from, y_from, z_from = mp_from.x0, mp_from.y0, mp_from.z0

        v_s_from = fun_s(x_from).reshape(-1, 1)

        v_v_from = fun_v(x_from, y_from, z_from)


        interface_from.set_variable_data(mp_name_from, var_s, v_s_from)
        interface_from.set_variable_data(mp_name_from, var_v, v_v_from)

        parameters_to = [{'model_part': mp_name_to, 'variables': [var_s, var_v]}]
        interface_to = data_structure.Interface(parameters_to, model)

        # check mapped values for 1D variable
        mapper((interface_from, mp_name_from, var_s),
               (interface_to, mp_name_to, var_s))
        v_s_to_ref = np.zeros((self.mp_to.size, 1))
        for i in range(self.mp_to.size):
            if self.mp_to.x0[i] < self.v_min or self.mp_to.x0[i] >self.v_max:
                v_s_to_ref[i] = 0
            else:
                v_s_to_ref[i] = v_s_from[i]

        v_s_to = interface_to.get_variable_data(mp_name_to, var_s)

        print("v_s_to")
        print(v_s_to)
        print("v_s_to_ref")
        print(v_s_to_ref)
        np.testing.assert_allclose(v_s_to, v_s_to_ref, rtol=1e-14)
        # check mapped values for 3D variable
        mapper((interface_from, mp_name_from, var_v),
               (interface_to, mp_name_to, var_v))
        v_v_to_ref = fun_v(x_to, y_to, z_to)
        for i in range(self.mp_to.size):
            if self.mp_to.x0[i] < self.v_min or self.mp_to.x0[i] >self.v_max:
                v_v_to_ref[i,0] = 0
                v_v_to_ref[i,1] = 0
            else:
                v_v_to_ref[i,0] = v_v_from[i,0]
                v_v_to_ref[i, 1] =v_v_from[i,1]
        v_v_to = interface_to.get_variable_data(mp_name_to, var_v)
        # print("v_v_to")
        # print(v_v_to)
        # print("v_v_to_ref")
        # print(v_v_to_ref)

        np.testing.assert_allclose(v_v_to, v_v_to_ref, rtol=1e-14)

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
                        pivot='tail', arrow_length_ratio=0.05, normalize=False, length=0.05, colors='r', linewidth=3)
            ax_v.quiver(x_to, y_to, z_to, v_v_to[:, 0], v_v_to[:, 1], v_v_to[:, 2],
                        pivot='tail', arrow_length_ratio=0.05, normalize=False, length=0.08)

            for ax in [ax_s, ax_v]:
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')

            plt.get_current_fig_manager().window.showMaximized()
            plt.show()
            plt.close()

if __name__ == '__main__':
    unittest.main()