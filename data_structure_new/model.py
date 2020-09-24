from coconut.data_structure_new.model_part import ModelPart


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

    """
    looping over ModelParts instead of ModelPart-names?
    --> requires new __iter__ and __next__
    --> when would we use __iter__ functionality??
        (mostly in Interface class?)
        
        leave it like this for now
    """

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
