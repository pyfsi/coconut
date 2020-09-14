# pyKratos imports
from .Element import Element
import math
import numpy as np

def Create(elem_id, nodes, gauss_nodes):
    return Quad4NodeElement(elem_id, nodes, gauss_nodes)


class Quad4NodeElement(Element):

    def __init__(self, elem_id, points, gauss_nodes):
        super(Quad4NodeElement, self).__init__(elem_id, points, gauss_nodes)
        g1 = 0.0
        g2 = 1 / math.sqrt(3)
        w1 = 2.0
        w2 = 1.0
        self.integration_local_points_weights_dict = {0:[],
                                                      1: [(np.array([g1, g1]), w1)],
                                                      4: [(np.array([g2, g2]), w2), (np.array([-1 * g2, g2]), w2),
                                                     (np.array([-1 * g2, -1 * g2]), w2), (np.array([g2, -1 * g2]), w2)]}

        if(len(self.GetPoints()) != 4):
            raise Exception("wrong number of nodes! should be 4!")

        # for point in self.GetPoints():
        #     if(point.Id < 0):
        #         raise Exception("point/node with Id smaller than 0 found")
        supported_number_gauss_nodes = self.integration_local_points_weights_dict.keys()
        nr_of_gauss_nodes = self.NumberOfGaussNodes()

        if nr_of_gauss_nodes not in self.integration_local_points_weights_dict.keys():
            raise Exception(f"{self.NumberOfGaussNodes()} gauss nodes not supported. Available options are: {supported_number_gauss_nodes}.")

        self.__gauss_node_id_int_wt_pair_dict = {}
        for gauss_node in self.GetGaussNodes():
            # if (gauss_node.Id < 0):  # *** doesn't work with string Id's
            #     raise Exception("gauss node with Id smaller than 0 found")
            index = self._FindGaussNodeIndex(gauss_node)
            integration_local_points_weights = self.integration_local_points_weights_dict[nr_of_gauss_nodes]
            self.__gauss_node_id_int_wt_pair_dict[gauss_node.Id] = integration_local_points_weights[index]



    def _FindGaussNodeIndex(self, node):
        nr_gauss_nodes =  self.NumberOfGaussNodes()
        integration_local_points_weights = self.integration_local_points_weights_dict[nr_gauss_nodes]
        for index, int_pt_wt_tuple in enumerate(integration_local_points_weights):
            local_gauss_point,_  = int_pt_wt_tuple
            global_gauss_point = self.GetGlobalCoordinate(local_gauss_point)
            distance = np.linalg.norm(global_gauss_point - node.Coordinates())
            if distance < 1e-5:  # TODO: change to relative criterion
                return index
        raise Exception(f"{node} is not a gauss point in the Quad4NodeElement element")




    def GetGlobalCoordinate(self,local_point):
        nodal_points = self.GetPoints()
        nodal_coordinates = np.array([point.Coordinates() for point in nodal_points])
        N = self._ShapeFunctions(local_point)
        global_coordinate = np.dot(N, nodal_coordinates)
        return global_coordinate



    def _ShapeFunctions(self, local_point):
        N = []
        N.append(0.25*(1.0 - local_point[0])*(1.0 - local_point[1]))
        N.append(0.25*(1.0 + local_point[0])*(1.0 - local_point[1]))
        N.append(0.25*(1.0 + local_point[0])*(1.0 + local_point[1]))
        N.append(0.25*(1.0 - local_point[0])*(1.0 + local_point[1]))

        return np.array(N)

    def _ShapeFunctionGradient(self, local_point):
        dN_dxi = [-0.25*1.0*(1.0 - local_point[1]), 0.25*1.0*(1.0 - local_point[1]), 0.25*1.0*(1.0 + local_point[1]), -0.25*1.0*(1.0 + local_point[1])]
        dN_eta = [-0.25*1.0*(1.0 - local_point[0]), -0.25*1.0*(1.0 + local_point[0]), 0.25*1.0*(1.0 + local_point[0]), 0.25*1.0*(1.0 - local_point[0])]
        shape_function_gradient =  [dN_dxi, dN_eta]
        return np.array(shape_function_gradient)

    def JacobianAtGaussNode(self, gauss_node):
        nodal_points = self.GetPoints()
        nodal_coordinates = np.array([point.Coordinates() for point in nodal_points])

        try:
            int_point, weight = self.__gauss_node_id_int_wt_pair_dict[gauss_node.Id]
        except KeyError:
            raise RuntimeError(f"{gauss_node} not found in the element")
        J = np.dot(self._ShapeFunctionGradient(int_point), nodal_coordinates)
        return J

    def _Jacobian(self, local_point):
        nodal_points = self.GetPoints()
        nodal_coordinates = np.array([point.Coordinates() for point in nodal_points])
        J_trans = np.dot(self._ShapeFunctionGradient(local_point), nodal_coordinates)
        return np.transpose(J_trans)


    def _DeterminantOfJacobian(self, local_point):
        J = self._Jacobian(local_point)
        det_J = math.sqrt(np.linalg.det((np.dot(np.transpose(J),J))))
        return det_J

    def Area(self):
        area = 0.0
        for gauss_node in self.GetGaussNodes():
            local_point, weight = self.__gauss_node_id_int_wt_pair_dict[gauss_node.Id]
            det_j = self._DeterminantOfJacobian(local_point)
            area += det_j*weight
        return area

    def TrueArea(self):
        area = 0.0
        integration_local_points_weights = self.integration_local_points_weights_dict[4]
        for int_point, weight in integration_local_points_weights:
            det_j = self._DeterminantOfJacobian(int_point)
            area += det_j * weight
        return area

    def CalculateElementAverage(self, variable):
        if variable.Type() == "Double":
            value = 0.0
            for gauss_node in self.GetGaussNodes():
                if gauss_node.SolutionStepsDataHas(variable):
                    local_point, weight = self.__gauss_node_id_int_wt_pair_dict[gauss_node.Id]
                    det_j = self._DeterminantOfJacobian(local_point)
                    value += det_j * weight * gauss_node.GetSolutionStepValue(variable)
            value /= self.Area()
        elif variable.Type() == "Array":
            value = np.zeros(3)
            for gauss_node in self.GetGaussNodes():
                if gauss_node.SolutionStepsDataHas(variable):
                    local_point, weight = self.__gauss_node_id_int_wt_pair_dict[gauss_node.Id]
                    det_j = self._DeterminantOfJacobian(local_point)
                    value += det_j * weight * np.array(gauss_node.GetSolutionStepValue(variable))
            value /= self.Area()
        else:
            raise RuntimeError(f"No implementation for {variable}")
        return value



