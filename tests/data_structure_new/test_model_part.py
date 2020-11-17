from coconut.data_structure_new.model_part import ModelPart
import unittest
import numpy as np

##Why default value of ids is None

class TestModelPart(unittest.TestCase):

    def setUp(self):
        self.correct_size = 3
        self.incorrect_size = 2
        self.name = 'MP1'
        self.ids = np.arange(0, self.correct_size)
        self.x0 = np.random.rand(self.correct_size)
        self.y0 = np.random.rand(self.correct_size)
        self.z0 = np.random.rand(self.correct_size)



    def test_instantiation(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)

    def test_model_part_name(self):
        self.name = 1
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)


    def test_coordinate_array(self):
        self.z0 = np.random.rand(self.incorrect_size)
        self.assertRaises(ValueError, ModelPart, self.name, self.x0,self.y0, self.z0, self.ids)

    def test_ids(self):
        self.ids = np.arange(0, self.incorrect_size)
        self.assertRaises(ValueError, ModelPart, self.name, self.x0, self.y0, self.z0, self.ids)
        # TODO: check if ids are integer
        #check for duplicity of ids
        self.ids = np.full(shape=self.correct_size, fill_value=1)
        self.assertRaises(ValueError, ModelPart,self. name, self.x0, self.y0, self.z0, self.ids)


    def test_size(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)
        self.assertEqual(mp.size, self.correct_size, f"Should be equal to {self.correct_size}")

    def test_attribute_change(self):
        mp = ModelPart(self.name, self.x0, self.y0, self.z0, self.ids)
        with self.assertRaises(AttributeError):
            mp.x0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.y0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.z0 = np.random.rand(self.correct_size)
        with self.assertRaises(AttributeError):
            mp.name = 'MP2'
        with self.assertRaises(AttributeError):
            mp.id = np.arange(0, self.correct_size)



if __name__ == '__main__':
    unittest.main()

