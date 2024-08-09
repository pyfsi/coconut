import os
import unittest
import warnings

from coconut.examples.post_processing.post_processing import *
from numpy.testing import assert_allclose


class TestPostProcessing(unittest.TestCase):
    gui = True

    def setUp(self):
        dir_name = os.path.realpath(os.path.dirname(__file__))  # path to fluent directory
        self.case_path = os.path.join(dir_name, 'case_results.pickle')

    def tearDown(self):
        for file in ('test.gif', 'test.png'):
            if file in os.listdir(os.curdir):
                os.remove(file)

    def test_data_selection(self):
        pp = PostProcess(self.case_path)
        print(pp)
        pp.print_summary()
        pp.print_info()
        with self.assertRaises(ValueError):
            pp.add_subset()
        with self.assertRaises(ValueError):
            pp.add_subset(model_part='wall')
        sx = pp.add_subset(interface='interface_x')
        sy = pp.add_subset(interface='interface_y')
        self.assertEqual(pp.get_subsets(), [sx, sy])

        self.assertEqual(sx.get_available_variables(), ['displacement', 'coordinates'])
        self.assertEqual(sy.get_available_variables(), ['pressure', 'traction'])
        assert_allclose(sx.get_all_times(), np.array([3, 4, 5, 6]) * 0.0001)
        assert_allclose(sx.get_all_steps(), np.array([3, 4, 5, 6]))
        initial_coordinates = sx.get_all_initial_coordinates()
        assert_allclose(initial_coordinates, np.array([
            [0., 0.005, -0.0225],
            [0., 0.005, -0.0175],
            [0., 0.005, -0.0125],
            [0., 0.005, -0.0075],
            [0., 0.005, -0.0025],
            [0., 0.005, 0.0025],
            [0., 0.005, 0.0075],
            [0., 0.005, 0.0125],
            [0., 0.005, 0.0175],
            [0., 0.005, 0.0225]]))
        mask1 = abs(initial_coordinates[:, 2]) < 0.005
        sx.select_points(mask1)
        mask2 = sx.get_all_times() <= 0.0004
        sx.select_times(mask2)

        initial_coordinates = sy.get_all_initial_coordinates()
        mask1 = (initial_coordinates[:, 2] > 0.01)
        sy.select_points(mask1)

        print(sx)
        print(sy)

        self.assertEqual(sx.get_size(), 2)
        self.assertEqual(sx.get_num_steps(), 1)
        assert_allclose(sx.get_initial_coordinates_selection(), np.array([[0., 0.005, -0.0025],
                                                                          [0., 0.005, 0.0025]]))
        assert_allclose(sx.get_times_selection(), np.array([0.0003, 0.0004]))
        assert_allclose(sx.get_steps_selection(), np.array([3, 4]))

        disp = np.array([[[48, 52, 56],
                          [60, 64, 68]],
                         [[49, 53, 57],
                          [61, 65, 69]]])
        assert_allclose(sx.get_values('displacement'), disp)
        assert_allclose(sx.get_values('displacement', component='x'), np.array(disp[:, :, 0]))
        assert_allclose(sx.get_values(), disp)

        pres = np.array([[[28], [32], [36]],
                      [[29], [33], [37]],
                      [[30], [34], [38]],
                      [[31], [35], [39]]])
        with self.assertRaises(ValueError):
            sy.get_values()
        assert_allclose(sy.get_values('pressure'), pres)

        sx.reset_times_selection()
        self.assertEqual(sx.get_num_steps(), 3)
        sx.reset_points_selection()
        self.assertEqual(sx.get_size(), 10)

    def test_visualization(self):
        pp = PostProcess(self.case_path)

        sx = pp.add_subset(interface='interface_x')
        sy = pp.add_subset(interface='interface_y')

        Animation2dDisplacement(sx, x_component='z', y_component='y').save('test.gif')
        Plot2dPressure(sy, x_component='z', time_step=3).save('test.png')
        Animation3dDisplacement(pp, component='y').save('test.gif')
        displacement_animation = Animation2dDisplacement(sx, x_component='z', y_component='y',
                                                         func_animation_settings=dict(frames=4, skip=1))
        print(displacement_animation)
        displacement_animation.save('test.gif')

        pressure_animation = Animation3dPressure(pp)
        pressure_animation.get_ax().set_title('Pressure animation')
        pressure_animation.get_figure().tight_layout()
        writer = ani.PillowWriter(fps=15, metadata=dict(artist='CoCoNuT'), bitrate=1800)
        pressure_animation.set_writer(writer)
        pressure_animation.save('test.gif')

        sx2 = pp.add_subset(interface='interface_x')
        coordinates_animation = Animation2dCoordinates(sx2, x_component='z', y_component='y', name='yz')
        coordinates_animation.remove_subset(0)

        line = coordinates_animation.get_line('yz')
        line.set(color='green', linestyle='', marker='o', markersize=5)
        coordinates_animation.save('test.gif')

        sx3 = pp.add_subset(interface='interface_x')
        sx3.select_points(sx3.get_all_initial_coordinates()[:, 2] > 0.02)
        time_evo = TimeEvolution(sx3, 'displacement')
        time_evo.add_subset(sx3, 'displacement', component='x')
        time_evo.save('test.png')

        if self.gui:
            plt.show()
            plt.close()


if __name__ == '__main__':
    unittest.main()
