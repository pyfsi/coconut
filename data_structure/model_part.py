import numpy as np


class ModelPart:

    def __init__(self, name, x0, y0, z0, id):
        """
        - inputs should have same size and be 1D
        - ModelPart contains no data, only coordinates
        - ModelPart has no idea: we use index in array
        """
        if type(name) is not str:
            raise ValueError('name should be a string')
        for coordinate in (x0, y0, z0):
            if type(coordinate) is not np.ndarray or coordinate.ndim != 1:
                raise ValueError('all coordinate arrays should be a 1D ndarray')
        if x0.size != y0.size or x0.size != z0.size:
            raise ValueError(f'all coordinates should have the same size:'
                             f'\n\tx0\t{x0.size}\n\ty0\t{y0.size}\n\tz0\t{z0.size}')
        if type(id) is not np.ndarray or id.dtype.kind not in np.typecodes['AllInteger'] or id.shape != x0.shape:
            raise ValueError(f'id array should be 1D array of integers of the same size as the coordinate arrays')
        _, counts = np.unique(id, return_counts=True)
        if np.any(counts > 1):
            raise ValueError(f'id array should be 1D array should have unique integer values')

        self.__name = name
        self.__x0 = x0.copy()
        self.__y0 = y0.copy()
        self.__z0 = z0.copy()
        self.__id = id.copy()

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
    def id(self):
        return self.__id.copy()

    @property
    def name(self):
        return self.__name

    @property
    def size(self):  # name is like numpy; or rename to number_of_points?
        return self.__x0.size

    def __eq__(self, other):
        if type(other) is ModelPart:
            return (self.__name == other.name and np.all(self.__x0 == other.x0) and np.all(self.__y0 == other.y0)
                    and np.all(self.__z0 == other.z0) and np.all(self.__id == other.id))
        return NotImplemented
