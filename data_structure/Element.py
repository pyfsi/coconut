from .data_value_container import DataValueContainer
from collections import OrderedDict

class Element(DataValueContainer):


    def __init__(self, elem_id, points, gauss_nodes):
        super(Element, self).__init__()
        self.Id = elem_id
        self.__points = points
        self.__gauss_nodes = gauss_nodes
        self.__variables = {}


    def GetPoint(self, point_index):
        return self.__points[point_index]

    def GetPoints(self):
        return self.__points

    def NumberOfPoints(self):
        return len(self.__points)


    def NumberOfGaussNodes(self):
        return len(self.__gauss_nodes)

    def GetGaussNodes(self):
        return self.__gauss_nodes


    def Initialize(self):
        pass

    def __str__(self):
        return  "Element #{} with {}".format(self.Id, self.__variables)

