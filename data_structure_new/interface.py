from coconut.data_structure_new.variables import variables_dimensions

import numpy as np
import copy


class Interface:

    def __init__(self, parameters, model):
        self.__model = model
        self.__parameters = parameters
        self.__data = {}

        # if not type(model) == Model:
        #     raise TypeError
        # if not type(parameters) == dict:
        #     raise TypeError

        for model_part_name in parameters:
            model_part = model.get_model_part(model_part_name)
            self.__data[model_part_name] = {}
            for variable in parameters[model_part_name]:
                if variable not in variables_dimensions:
                    raise ValueError(f'Invalid variable name "{variable}"')
                shape = (model_part.size, variables_dimensions[variable])
                self.__data[model_part_name][variable] = np.zeros(shape)

        tmp = []
        for model_part_name in self.parameters:
            for variable in self.parameters[model_part_name]:
                tmp.append((model_part_name, variable))
        self.__model_part_variable_pairs = tmp

    @property
    def model_part_variable_pairs(self):
        return copy.deepcopy(self.__model_part_variable_pairs)

    @property
    def parameters(self):
        return copy.deepcopy(self.__parameters)

    @property
    def size(self):
        s = 0
        for model_part_name, variable in self.model_part_variable_pairs:
            s += self.__data[model_part_name][variable].size
        return s

    def copy_old(self):
        # create new Interface
        interface = Interface(self.__parameters, self.__model)

        # copy data
        interface.set_interface_data(self.get_interface_data())

        return interface

    def copy(self):
        # create new Interface
        interface = Interface(self.__parameters, self.__model)

        # copy data
        interface += self  # uses fast __iadd__ method to transfer data

        return interface

    def __repr__(self):
        repr = 'Interface that refers to ModelParts'
        for model_part_name in self.parameters:
            model_part = self.__model.get_model_part(model_part_name)
            repr += f'\n\t"{model_part.name}" of size {model_part.size} with variables'
            for variable in self.__data[model_part_name]:
                repr += f'\n\t\t{variable} with {variables_dimensions[variable]} components'
        return repr

    def get_variable_data(self, model_part_name, variable):
        # *** always returns copies!
        if (model_part_name, variable) not in self.model_part_variable_pairs:
            raise KeyError
        return self.__data[model_part_name][variable].copy()

    def set_variable_data(self, model_part_name, variable, data):
        # *** this changes the original data!
        if not isinstance(data, np.ndarray):
            raise ValueError(f'data is of type {type(data)}, but must be ndarray')
        shape = self.__data[model_part_name][variable].shape
        if data.shape != shape:
            raise ValueError(f'ndarray has shape {data.shape} instead of shape {shape}')
        self.__data[model_part_name][variable] = data.copy()

    def get_interface_data(self):
        data = np.empty(0)
        for model_part_name, variable in self.model_part_variable_pairs:
            data = np.concatenate((data, self.get_variable_data(model_part_name, variable).flatten()))
        return data

    def set_interface_data(self, data):
        if not isinstance(data, np.ndarray):
            raise ValueError(f'data is of type {type(data)}, but must be ndarray')
        if data.shape != (self.size,):
            raise ValueError(f'ndarray has shape {data.shape} instead of shape {(self.size,)}')
        index = 0
        for model_part_name, variable in self.model_part_variable_pairs:
            tmp = self.get_variable_data(model_part_name, variable)
            self.set_variable_data(model_part_name, variable,
                                   data[index:index + tmp.size].reshape(tmp.shape))
            index += tmp.size

    """
    For one type of operation (add), I implemented the 3 required methods.
    They do:
        __add__:  self + other
        __radd__: other + self
        __iadd__: self += other
    I also implemented 3 different ways to do these operations (from high level 
    operations to low level operations). 
    The time I measured is given in the methods. 
    Conclusions:
        - low-level methods are much faster
        - "+=" is much faster than "+"
        - copy() is much faster if "+=" is used instead of get/set_interface_data
    The low-level one is faster, as expected.
    """

    def __add__(self, other):
        """
        speed for case
            1: 375
            2: 270
            3: 189
        with new copy() method:
            1: 256
            2: 152
            3: 69
        """
        case = 3
        if case == 1:
            result = self.copy()
            if isinstance(other, Interface):
                # loop over 1 Interface.model_part_variable_pairs --> gives automatic check if Interface are same
                result.set_interface_data(self.get_interface_data() + other.get_interface_data())
            elif isinstance(other, (int, float, np.integer, np.floating)):
                result.set_interface_data(self.get_interface_data() + other)
            else:
                return NotImplemented
            return result
        if case == 2:
            result = self.copy()
            if isinstance(other, Interface):
                for model_part_name, variable in self.model_part_variable_pairs:
                    result.set_variable_data(model_part_name, variable,
                                             self.get_variable_data(model_part_name, variable) +
                                             other.get_variable_data(model_part_name, variable))
            else:
                return NotImplemented
            return result
        if case == 3:
            result = self.copy()
            if isinstance(other, Interface):
                for model_part_name, variable in self.model_part_variable_pairs:
                    result.__data[model_part_name][variable] += self.__data[model_part_name][variable]
            else:
                return NotImplemented
            return result

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):  # iadd does not need to use copy
        """
        speed for case
            1: 107
            2: 35
            3: 18
        """
        case = 3
        if isinstance(other, Interface):
            if case == 1:
                self.set_interface_data(self.get_interface_data() + other.get_interface_data())
            if case == 2:
                for model_part_name, variable in self.model_part_variable_pairs:
                    self.set_variable_data(model_part_name, variable,
                                           self.get_variable_data(model_part_name, variable) +
                                           other.get_variable_data(model_part_name, variable))
            if case == 3:
                for model_part_name, variable in self.model_part_variable_pairs:
                    self.__data[model_part_name][variable] += other.__data[model_part_name][variable]
        else:
            return NotImplemented
        return self
