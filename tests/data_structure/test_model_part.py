from coconut.data_structure.model_part import ModelPart

import unittest
import numpy as np


class TestModelPart(unittest.TestCase):

    def setUp(self):
        self.correct_size = 3
        self.incorrect_size = 2
        self.name = 'mp1'
        self.ids = np.arange(self.correct_size)
        np.random.shuffle(self.ids)
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

    def test_eq(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)
        name_2 = 'mp2'
        x0_2 = self.x0.copy()
        x0_2[np.random.randint(self.correct_size)] = np.random.rand()
        y0_2 = self.y0.copy()
        y0_2[np.random.randint(self.correct_size)] = np.random.rand()
        z0_2 = self.z0.copy()
        z0_2[np.random.randint(self.correct_size)] = np.random.rand()
        ids_2 = self.ids.copy()
        ids_2[np.random.randint(self.correct_size)] = np.max(self.ids) + np.random.randint(1, 100)

        args = (self.name, self.x0, self.y0, self.z0, self.ids)
        mp2 = ModelPart(*args)
        self.assertEqual(mp, mp2)

        changed_arguments = [name_2, x0_2, y0_2, z0_2, ids_2]
        for i in range(5):
            args_2 = list(args)
            args_2[i] = changed_arguments[i]
            mp2 = ModelPart(*tuple(args_2))
            self.assertNotEqual(mp, mp2)


if __name__ == '__main__':
    unittest.main()
