class ModelPart:

    def __init__(self, name, x0, y0, z0):
        """
        - inputs should have same size and be 1D
        - ModelPart contains no data, only coordinates
        - ModelPart has no idea: we use index in array
        """
        # TODO: check sizes and shapes of matrices, raise error
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
    def size(self):  # name is like numpy; or rename to number_of_points?
        return self.__x0.size

