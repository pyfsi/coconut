from coconut.data_structure.model import Model
from coconut.data_structure.model_part import ModelPart
from coconut.data_structure.interface import Interface
from coconut.coupling_components.tools import create_instance

import numpy as np


# create test-stuff for interface mapper
model_part_size = 3
parameters = {'interface_from':
                  [{'model_part': 'mp1', 'variables': ['pressure', 'traction']},
                   {'model_part': 'mp2', 'variables': ['displacement']}],
              'interface_to':
                  [{'model_part': 'mp3', 'variables': ['pressure', 'traction']},
                   {'model_part': 'mp4', 'variables': ['displacement']}]}
model = Model()

n = [2, 3, 4, 5]
for i in range(4):
    x0 = np.linspace(0, 1, n[i])
    y0, z0 = np.zeros(n[i]), np.zeros(n[i])
    model.create_model_part(f'mp{i + 1}', x0, y0, z0, np.arange(n[i]))

interface_from = Interface(parameters['interface_from'], model)
interface_to = Interface(parameters['interface_to'], model)

print(interface_from)
print(interface_to)

# create Interface mapper
tmp = {'type': 'mappers.nearest', 'settings': {'directions': ['x', 'y']}}
mapper_parameters = {'type': 'mappers.interface', 'settings': tmp}

mapper = create_instance(mapper_parameters)
mapper.Initialize(interface_from, interface_to)

mapper(interface_from, interface_to)

mapper.PrintInfo('')

print('------------------')


# create test-stuff for nearest neighbour mapper
model = Model()

n = [7, 11]
mps = []
for i in range(2):
    x0 = np.linspace(0, 1, n[i])
    y0, z0 = np.zeros(n[i]), np.zeros(n[i])
    mps.append(ModelPart(f'mp{i + 1}', x0, y0, z0, np.arange(n[i])))
mp1, mp2 = mps

# create ModelPart mapper
mapper_parameters = {"type": "mappers.nearest",
                     "settings": {"directions": ["x", "y"]}}

mapper = create_instance(mapper_parameters)
mapper.Initialize(mp1, mp2)





