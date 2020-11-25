from coconut.data_structure_new.model_part import ModelPart

import unittest
import numpy as np


class TestModelPart(unittest.TestCase):

    def setUp(self):
        self.correct_size = 3
        self.incorrect_size = 2
        self.name = 'mp1'
        self.ids = np.arange(self.correct_size)
        self.x0 = np.random.rand(self.correct_size)
        self.y0 = np.random.rand(self.correct_size)
        self.z0 = np.random.rand(self.correct_size)

    def test_instantiation(self):
        ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)

    def test_model_part_name(self):
        self.name = 1
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)

    def test_coordinate_array(self):
        self.z0 = np.random.rand(self.incorrect_size)
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)

    def test_ids(self):
        # check ids have correct size
        self.ids = np.arange(self.incorrect_size)
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)
        # check ids are integers
        self.ids = np.arange(self.correct_size) + 0.5
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)
        # check id is 1d numpy array
        self.ids = list(range(self.correct_size))
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)
        self.ids = np.arange(self.correct_size).reshape(1, -1)
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)
        # check for duplicity of ids
        self.ids = np.full(shape=self.correct_size, fill_value=1)
        self.assertRaises(ValueError, ModelPart, self. name, self.x0, self.y0, self.z0, self.ids)

    def test_size(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)
        self.assertEqual(mp.size, self.correct_size, f"should be equal to {self.correct_size}")

    def test_attribute_change(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)
        with self.assertRaises(AttributeError):
            mp.x0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.y0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.z0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.name = 'mp2'
        with self.assertRaises(AttributeError):
            mp.id = np.arange(self.correct_size)


if __name__ == '__main__':
    unittest.main()
