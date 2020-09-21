from coconut.coupling_components.tools import quicktimer
from coconut.data_structure.Model import Model as ModelOld
from coconut.coupling_components.interface import Interface as InterfaceOld
from coconut.data_structure.Parameters import Parameters
from coconut import data_structure

import numpy as np
import copy
import json


# >> new code

class Model:
    def __init__(self):
        self.__model_parts = {}

    def create_model_part(self, name, x0, y0, z0):
        if name in self.__model_parts:
            raise ValueError  # TODO
        self.__model_parts[name] = ModelPart(name, x0, y0, z0)

    def get_model_part(self, name):
        if name not in self.__model_parts:
            raise ValueError  # TODO
        return self.__model_parts[name]

    def __iter__(self):  # iterate over names of ModelParts
        return iter(self.__model_parts)

    def __next__(self):
        return next(self.__model_parts)

    def __repr__(self):
        repr = 'Model that consists of ModelParts'
        for model_part_name in self.__model_parts:
            repr += (f'\n\t"{self.__model_parts[model_part_name].name}" ' +
                    f'of size {self.__model_parts[model_part_name].size}')
        return repr


class ModelPart:

    def __init__(self, name, x0, y0, z0):
        """
        - inputs should have same size and be 1D
        - ModelPart contains no data, only coordinates
        - ModelPart has no idea: we use index in array
        """
        self.__name = name
        self.__x0 = x0.copy()
        self.__y0 = y0.copy()
        self.__z0 = z0.copy()

    def __repr__(self):
        return f'ModelPart "{self.name}" of size {self.size}'

    @property
    def x0(self):
        return self.__x0.copy()

    @property
    def y0(self):
        return self.__y0.copy()

    @property
    def z0(self):
        return self.__z0.copy()

    @property
    def name(self):
        return self.__name

    @property
    def size(self):
        return self.__x0.size


class Interface:

    def __init__(self, parameters=None, model=None):

        self.__model_parts = {}
        self.__data = {}

        if parameters and model:
            self.__parameters = parameters
            for model_part_name in parameters:
                model_part = model.get_model_part(model_part_name)
                self.__model_parts[model_part_name] = model_part
                self.__data[model_part_name] = {}
                for variable in parameters[model_part_name]:
                    if variable not in variables_components:
                        raise ValueError(f'Invalid variable name "{variable}"')
                    shape = (model_part.size, variables_components[variable])
                    self.__data[model_part_name][variable] = np.zeros(shape)

    def key_pairs(self):
        # to shorten double loops over model_part_names and variables
        out = []
        for model_part_name in self.parameters:
            for variable in self.parameters[model_part_name]:
                out.append((model_part_name, variable))
        return tuple(out)

    @property
    def parameters(self):
        return copy.deepcopy(self.__parameters)

    @property
    def size(self):
        s = 0
        for model_part_name, variable in self.key_pairs():
            s += self.__data[model_part_name][variable].size
        return s

    def copy(self):
        # create new Interface
        interface = Interface()

        # reference ModelParts
        interface._Interface__parameters = self.__parameters
        interface._Interface__model_parts = self.__model_parts

        # copy data
        for model_part_name, variable in self.key_pairs():
            if model_part_name not in interface._Interface__data:
                interface._Interface__data[model_part_name] = {}
            interface._Interface__data[model_part_name][variable] = \
                self.__data[model_part_name][variable].copy()

        return interface

    def __repr__(self):
        repr = 'Interface that refers to ModelParts'
        for model_part_name in self.__model_parts:
            model_part = self.__model_parts[model_part_name]
            repr += f'\n\t"{model_part.name}" of size {model_part.size} with variables'
            for variable in self.__data[model_part_name]:
                repr += f'\n\t\t{variable} with {variables_components[variable]} components'
        return repr

    def get_variable_data(self, model_part_name, variable):
        # always returns copies!
        # TODO: check if model_part and variable exist?
        return self.__data[model_part_name][variable].copy()

    def set_variable_data(self, model_part_name, variable, data):
        # this changes the original data!
        if not isinstance(data, np.ndarray):
            raise ValueError(f'data is of type {type(data)}, but must be ndarray')
        shape = self.__data[model_part_name][variable].shape
        if data.shape != shape:
            raise ValueError(f'ndarray has shape {data.shape} instead of shape {shape}')
        self.__data[model_part_name][variable] = data.copy()

    def get_interface_data(self):
        # like getNumpyArray
        lst = []
        for model_part_name, variable in self.key_pairs():
            lst.append(self.get_variable_data(model_part_name, variable).flatten())
        return np.concatenate(tuple(lst))

    def set_interface_data(self, data):
        # like setNumpyArray
        if not isinstance(data, np.ndarray):
            raise ValueError(f'data is of type {type(data)}, but must be ndarray')
        if data.shape != (self.size,):
            raise ValueError(f'ndarray has shape {data.shape} instead of shape {(self.size,)}')
        index = 0
        for model_part_name, variable in self.key_pairs():
            tmp = self.get_variable_data(model_part_name, variable)
            self.set_variable_data(model_part_name, variable,
                                   data[index:index + tmp.size].reshape(tmp.shape))
            index += tmp.size

    def __add__(self, other):
        # faster implementation possible I guess...
        result = self.copy()
        if isinstance(other, Interface):
            result.set_interface_data(self.get_interface_data() + other.get_interface_data())
        elif isinstance(other, (int, float, np.integer, np.floating)):
            result.set_interface_data(self.get_interface_data() + other)
        else:
            return NotImplemented
        return result


variables_components = {}
variables_components['pressure'] = 1
variables_components['traction'] = 3



# >> testing code


model = Model()
for i in range(3):
    coords = np.random.rand(5 + i, 3)
    model.create_model_part(f'mp{i}', coords[:, 0], coords[:, 1], coords[:, 2])

print(model)
mp1 = model.get_model_part('mp1')

print(mp1)

par = {'mp0': ['pressure'],
       'mp1': ['pressure', 'traction']}  # *** order correct?? orderd dicts??

interface = Interface(par, model)
print(interface)

print(interface.get_interface_data())

interface.set_variable_data('mp1', 'pressure', np.random.rand(6, 1))

print(interface.get_variable_data('mp1', 'pressure'))
print(interface.get_interface_data())
print('\n' * 5)



# >> test speed
n = 10000
coords = np.random.rand(n, 3)

par = {'mp0': ['pressure'],
       'mp1': ['pressure', 'traction']}
tmp = {'mp0': ['PRESSURE'],
       'mp1': ['PRESSURE', 'TRACTION']}
par_old = Parameters(json.dumps(tmp))

# setup with new data structure
model = Model()
for model_part_name in par:
    model.create_model_part(model_part_name, *np.hsplit(coords, 3))
interface = Interface(par, model)

# setup with old data structure
model_old = ModelOld()
for model_part_name in par_old.keys():
    mp = model_old.CreateModelPart(model_part_name)
    for var in par_old[model_part_name].list():
        mp.AddNodalSolutionStepVariable(vars(data_structure)[var.GetString()])
    for i in range(n):
        mp.CreateNewNode(i, coords[i, 0], coords[i, 1], coords[i, 2])
interface_old = InterfaceOld(model_old, par_old)

# comparison
print('\nget and set numpy array in interface')
with quicktimer('new', t=1, ms=True):
    data = interface.get_interface_data()
    interface.set_interface_data(np.random.rand(*data.shape))
with quicktimer('old', t=1, ms=True):
    data_old = interface_old.GetNumpyArray()
    interface_old.SetNumpyArray(np.random.rand(*data_old.shape))

print('\ncopy interface')
with quicktimer('new', t=1, ms=True):
    interface_2 = interface.copy()
with quicktimer('old', t=1, ms=True):
    interface_old_2 = interface_old.deepcopy()

print('\nadd interfaces')
with quicktimer('new', t=1, ms=True):
    interface + interface_2
with quicktimer('old', t=1, ms=True):
    interface_old + interface_old_2
