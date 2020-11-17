from coconut.data_structure_new.model import Model
import unittest
import numpy as np

##Why default value of ids is None

class TestModel(unittest.TestCase):

    def setUp(self): # this function is called after every test
        self.model = Model()
        self.size = 3
        self.name = 'MP1'
        self.ids = np.arange(0, self.size)
        self.x0 = np.random.rand(self.size)
        self.y0 = np.random.rand(self.size)
        self.z0 = np.random.rand(self.size)
        self.model.create_model_part(self.name, self.x0, self.y0, self.z0, self.ids)



    def test_create_model_part(self):

        name = 'MP1'
        ids = np.arange(self.size-1, 2*self.size)
        x0 = np.random.rand(self.size)
        y0 = np.random.rand(self.size)
        z0 = np.random.rand(self.size)
        self.assertRaises(ValueError, self.model.create_model_part, name, x0, y0, z0, ids)

    def test_get_model_part(self):
        self.model.get_model_part('MP1') # MP1 is already present (from setUp method)
        self.assertRaises(ValueError, self.model.get_model_part, 'MP2')

    def test_model_part_iterator(self):
        name_list = [self.name]
        name = 'MP2'
        name_list.append(name)
        ids = np.arange(0, self.size)
        x0 = np.random.rand(self.size)
        y0 = np.random.rand(self.size)
        z0 = np.random.rand(self.size)
        self.model.create_model_part(name, x0, y0, z0, ids)

        for i, model_part_name in enumerate(self.model):
            self.assertEqual(model_part_name, name_list[i])



if __name__ == '__main__':
    unittest.main()