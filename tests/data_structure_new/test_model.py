from coconut.data_structure_new.model import Model
from coconut.data_structure_new.model_part import ModelPart

import unittest
import numpy as np


class TestModel(unittest.TestCase):

    def setUp(self):  # this function is called before every test
        self.model = Model()
        self.size = 3
        self.name = 'mp1'
        self.ids = np.arange(self.size)
        self.x0 = np.random.rand(self.size)
        self.y0 = np.random.rand(self.size)
        self.z0 = np.random.rand(self.size)
        self.model.create_model_part(self.name, self.x0, self.y0, self.z0, self.ids)

    def test_create_model_part(self):
        name = 'mp1'
        ids = np.arange(self.size)
        x0 = np.random.rand(self.size)
        y0 = np.random.rand(self.size)
        z0 = np.random.rand(self.size)
        self.assertRaises(ValueError, self.model.create_model_part, name, x0, y0, z0, ids)
        name = 'mp2'
        size = 5
        ids = np.arange(size)
        x0 = np.random.rand(size)
        y0 = np.random.rand(size)
        z0 = np.random.rand(size)
        self.assertIsInstance(self.model.create_model_part(name, x0, y0, z0, ids), ModelPart)

    def test_get_model_part(self):
        self.model.get_model_part('mp1')  # mp1 is already present (from setUp method)
        self.assertRaises(ValueError, self.model.get_model_part, 'mp2')

    def test_model_part_iterator(self):
        name_list = [self.name]
        name = 'mp2'
        name_list.append(name)
        ids = np.arange(self.size)
        x0 = np.random.rand(self.size)
        y0 = np.random.rand(self.size)
        z0 = np.random.rand(self.size)
        self.model.create_model_part(name, x0, y0, z0, ids)

        for i, model_part_name in enumerate(self.model):
            self.assertEqual(model_part_name, name_list[i])

    def test_attribute_change(self):
        with self.assertRaises(AttributeError):
            self.model.model_parts


if __name__ == '__main__':
    unittest.main()
