from coconut.data_structure.model_part import ModelPart


class Model:
    def __init__(self):
        self.__model_parts = {}

    def create_model_part(self, name, x0, y0, z0, id):
        if name in self.__model_parts:
            raise ValueError(f'model already has model part with name "{name}"')
        self.__model_parts[name] = ModelPart(name, x0, y0, z0, id)
        return self.__model_parts[name]

    def get_model_part(self, name):
        if name not in self.__model_parts:
            raise ValueError(f'no model part with name "{name}" in model')
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

    def __eq__(self, other):
        if type(other) is Model:
            return self.__model_parts == other.__model_parts
        return NotImplemented
