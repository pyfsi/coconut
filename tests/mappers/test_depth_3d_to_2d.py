from coconut import data_structure
from coconut.tools import create_instance
from coconut.tests.mappers import test_depth_2d_to_3d

import numpy as np
import unittest


class TestMapperDepth3DTo2D(test_depth_2d_to_3d.TestMapperDepth2DTo3D):
    gui = False

    def setUp(self):
        self.parameters = {'type': 'mappers.depth_3d_to_2d',
                           'settings':
                               {'direction_depth': 'z',
                                'coordinates_depth': [-1, 1, -2, 2]}
                           }
        self.forward = False

    def test_call(self):
        def fun_s(x):
            return 1. + 2.5 * x

        def fun_v(x, y, z):
            v_x = 0.1 + 2.5 * x
            v_y = 0.1 + 2.5 * y
            v_z = -4 * z
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
            # import here to avoid error on systems without matplotlib
            import matplotlib.pyplot as plt
            from matplotlib import cm

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


if __name__ == '__main__':
    unittest.main()
