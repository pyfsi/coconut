from coconut import data_structure
from coconut.coupling_components.tools import check_bounding_box

import unittest
import numpy as np
import warnings


def split(coords):
    x, y, z = np.hsplit(coords, 3)
    return x.flatten(), y.flatten(), z.flatten()


class TestMapperInterpolator(unittest.TestCase):
    def test_bounding_box_1d(self):
        model = data_structure.Model()

        with self.assertWarns(Warning) as tmp:  # check that no warnings are thrown
            coords = np.array([[0, 0, 0], [0, 0, 1]])
            mp_a = model.create_model_part('mp_a', *split(coords), np.arange(2))
            mp_b = model.create_model_part('mp_b_1', *split(coords), np.arange(2))
            check_bounding_box(mp_a, mp_b)

            coords = np.vstack((coords, np.array([[0, 0, 1.01]])))
            mp_b = model.create_model_part('mp_b_2', *split(coords), np.arange(3))
            check_bounding_box(mp_a, mp_b)

            coords = np.vstack((coords, np.array([[0, 0, -.01]])))
            mp_b = model.create_model_part('mp_b_3', *split(coords), np.arange(4))
            check_bounding_box(mp_a, mp_b)

            warnings.warn('msg', Warning)
        self.assertEqual(len(tmp.warnings), 1)

        coords = np.vstack((coords, np.array([[0, 0, 1.1]])))
        mp_b = model.create_model_part('mp_b_4', *split(coords), np.arange(5))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[0, 0, 1.25]])))
        mp_b = model.create_model_part('mp_b_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[0, 0, -.25]])))
        mp_b = model.create_model_part('mp_b_6', *split(coords), np.arange(7))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[0, 0, 2], [0, 0, -1]])))
        mp_b = model.create_model_part('mp_b_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

    def test_bounding_box_2d(self):
        model = data_structure.Model()

        with self.assertWarns(Warning) as tmp:  # check that no warnings are thrown
            coords = np.array([[0, 0, 0], [1, 0, 1]])
            mp_a = model.create_model_part('mp_a', *split(coords), np.arange(2))
            mp_b = model.create_model_part('mp_b_1', *split(coords), np.arange(2))
            check_bounding_box(mp_a, mp_b)
            
            coords = np.vstack((coords, np.array([[1.01, 0, 1.01]])))
            mp_b = model.create_model_part('mp_b_2', *split(coords), np.arange(3))
            check_bounding_box(mp_a, mp_b)
    
            coords = np.vstack((coords, np.array([[-.01, 0, -.01]])))
            mp_b = model.create_model_part('mp_b_3', *split(coords), np.arange(4))
            check_bounding_box(mp_a, mp_b)
            
            warnings.warn('msg', Warning)
        self.assertEqual(len(tmp.warnings), 1)

        coords = np.vstack((coords, np.array([[1.1, 0, 1.1]])))
        mp_b = model.create_model_part('mp_b_4', *split(coords), np.arange(5))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[1.25, 0, 1.25]])))
        mp_b = model.create_model_part('mp_b_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[-.25, 0, -.25]])))
        mp_b = model.create_model_part('mp_b_6', *split(coords), np.arange(7))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[2, 0, 2], [-1, 0, -1]])))
        mp_b = model.create_model_part('mp_b_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

    def test_bounding_box_3d(self):
        model = data_structure.Model()

        with self.assertWarns(Warning) as tmp:  # check that no warnings are thrown
            coords = np.array([[0, 0, 0], [1, 1, 1]])
            mp_a = model.create_model_part('mp_a', *split(coords), np.arange(2))
            mp_b = model.create_model_part('mp_b_1', *split(coords), np.arange(2))
            check_bounding_box(mp_a, mp_b)

            coords = np.vstack((coords, np.array([[1.01, 1.01, 1.01]])))
            mp_b = model.create_model_part('mp_b_2', *split(coords), np.arange(3))
            check_bounding_box(mp_a, mp_b)

            coords = np.vstack((coords, np.array([[-.01, -.01, -.01]])))
            mp_b = model.create_model_part('mp_b_3', *split(coords), np.arange(4))
            check_bounding_box(mp_a, mp_b)

            warnings.warn('msg', Warning)
        self.assertEqual(len(tmp.warnings), 1)

        coords = np.vstack((coords, np.array([[1.1, 1.1, 1.1]])))
        mp_b = model.create_model_part('mp_b_4', *split(coords), np.arange(5))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[1.25, 1.25, 1.25]])))
        mp_b = model.create_model_part('mp_b_5', *split(coords), np.arange(6))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[-.25, -.25, -.25]])))
        mp_b = model.create_model_part('mp_b_6', *split(coords), np.arange(7))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.vstack((coords, np.array([[2, 2, 2], [-1, -1, -1]])))
        mp_b = model.create_model_part('mp_b_7', *split(coords), np.arange(9))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

    def test_bounding_box_lines(self):
        # check if method works for lines aligned with coordinate axes in 2D
        model = data_structure.Model()

        with self.assertWarns(Warning) as tmp:  # check that no warnings are thrown
            coords = np.array([[0, 1, 0], [1, 1, 0]])
            mp_a = model.create_model_part('mp_a', *split(coords), np.arange(2))

            coords = np.array([[0, 1, 0], [1.01, 1.01, 0.]])
            mp_b = model.create_model_part('mp_b_1', *split(coords), np.arange(2))
            check_bounding_box(mp_a, mp_b)

            warnings.warn('msg', Warning)
        self.assertEqual(len(tmp.warnings), 1)

        coords = np.array([[0, 1.05, 0], [1, 1.05, 0]])
        mp_b = model.create_model_part('mp_b_2', *split(coords), np.arange(2))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.array([[0, 1.25, 0], [1, 1.25, 0]])
        mp_b = model.create_model_part('mp_b_3', *split(coords), np.arange(2))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

        coords = np.array([[0, .85, 0], [1, 1.15, 0]])
        mp_b = model.create_model_part('mp_b_4', *split(coords), np.arange(2))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

    def test_bounding_box_planes(self):
        # check if method works for planes aligned with coordinate axes in 3D
        model = data_structure.Model()

        with self.assertWarns(Warning) as tmp:  # check that no warnings are thrown
            coords = np.array([[0, 1, 0], [1, 1, 0]])
            mp_a = model.create_model_part('mp_a', *split(coords), np.arange(2))

            coords = np.array([[0, 1, 0], [1.01, 1.01, 0.]])
            mp_b = model.create_model_part('mp_b_1', *split(coords), np.arange(2))
            check_bounding_box(mp_a, mp_b)

            warnings.warn('msg', Warning)
        self.assertEqual(len(tmp.warnings), 1)

        coords = np.array([[0, 1.05, 0], [1, 1.05, 0]])
        mp_b = model.create_model_part('mp_b_2', *split(coords), np.arange(2))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))

        coords = np.array([[0, 1.25, 0], [1, 1.25, 0]])
        mp_b = model.create_model_part('mp_b_3', *split(coords), np.arange(2))
        self.assertRaises(ValueError, check_bounding_box, *(mp_a, mp_b))

        coords = np.array([[0, .85, 0], [1, 1.15, 0]])
        mp_b = model.create_model_part('mp_b_4', *split(coords), np.arange(2))
        self.assertWarns(Warning, check_bounding_box, *(mp_a, mp_b))


if __name__ == '__main__':
    unittest.main()
