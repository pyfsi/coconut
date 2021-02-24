from coconut.tools import create_instance
from coconut.coupling_components.component import Component
from coconut import tools
from coconut import data_structure

import unittest
import numpy as np
import json

class Test_mapped_update_tesimplementations(unittest.TestCase):
    def test_mapper_combined(self):
        with open('mappers/test_mapped_udate_testimplementations.json') as parameter_file:
            parameters = json.load(parameter_file)

        mp_name_fluid = 'fluid'
        mp_name_structure = 'structure'


        self.model_fluid = data_structure.Model()
        self.model_structure = data_structure.Model()

        n_in = 10

        n_from = n_in

        x_fluid = np.linspace(-0.5, 3, n_in)
        r_fluid = 1 + 0.07 * np.sin(x_fluid * 600)

        x_in = np.zeros(n_from)
        y_in = np.zeros(n_from)
        z_in = np.zeros(n_from)

        i = 0
        for k in range(n_in):
            x_in[i] = x_fluid[k]
            y_in[i] = r_fluid[k] * np.cos(np.radians(2.5))
            z_in[i] = 0
            i += 1

        self.model_fluid.create_model_part(mp_name_fluid, x_in, y_in, z_in, np.arange(n_from))
        self.parameters_input_from = [{'model_part':mp_name_fluid, 'variables':['pressure', 'traction']}]
        self.parameters_output_to = [{'model_part': mp_name_fluid, 'variables': ['displacement']}]
        interface_input_from = data_structure.Interface(self.parameters_input_from, self.model_fluid)
        interface_output_to = data_structure.Interface(self.parameters_output_to, self.model_fluid)

        x_structure = np.linspace(-0.5, 3, 15)
        r_structure = 1 + 0.07 * np.sin(x_fluid * 600)

        x_in = np.zeros(15)
        y_in = np.zeros(15)
        z_in = np.zeros(15)

        i = 0
        for k in range(n_in):
            x_in[i] = x_structure[k]
            y_in[i] = r_structure[k]
            z_in[i] = 0
            i += 1