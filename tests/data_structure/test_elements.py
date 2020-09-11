from coconut import data_structure
from coconut.data_structure import KratosUnittest
import numpy as np
import math

### TODO: implement 1d line element

class TestQuad4NodeElement(KratosUnittest.TestCase):

    def test_element_with_without_gauss_nodes(self):
        g1 = 0.0
        g2 = (1 / math.sqrt(3))
        model = data_structure.Model()
        model_part = model.CreateModelPart('test_model_part')
        model_part.AddNodalSolutionStepVariable(vars(data_structure)['PRESSURE'])
        model_part.CreateNewNode(1, -1.0, -1.0, 0.0)
        model_part.CreateNewNode(2, 1.0, -1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, -1.0, 1.0, 0.0)
        model_part.CreateNewNode(5, g2, g2, 0.0)
        model_part.CreateNewNode(6, -1 * g2, g2, 0.0)
        model_part.CreateNewNode(7, -1 * g2, -1 * g2, 0.0)
        model_part.CreateNewNode(8, g2, -1 * g2, 0.0)
        model_part.CreateNewNode(9, g1, g1, 0.0)
        elem_1 = model_part.CreateNewElement('Quad4NodeElement', 1, [1, 2, 3, 4], [5,6,7,8])
        elem_2 = model_part.CreateNewElement('Quad4NodeElement', 2, [1, 2, 3, 4], [9])
        elem_3 = model_part.CreateNewElement('Quad4NodeElement', 3, [1, 2, 3, 4])


        area_gauss_4 = elem_1.Area()
        area_gauss_1 = elem_2.Area()
        area_gauss_0 = elem_3.Area()

        self.assertEqualTolerance(area_gauss_4, 4.0,1e-10)
        self.assertEqualTolerance(area_gauss_1, 2.0, 1e-10)
        self.assertEqualTolerance(area_gauss_0, 0.0, 1e-10)

    def test_element_integration(self):
        model = data_structure.Model()
        v_1 = np.array([2, 1, 1])
        v_2 = np.array([1, 2, 0.5])
        v = v_1+v_2
        g = (v_1 + v_2)/2
        model_part = model.CreateModelPart('test_model_part')
        model_part.AddNodalSolutionStepVariable(vars(data_structure)['PRESSURE'])
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, v_1[0], v_1[1], v_1[2])
        model_part.CreateNewNode(3, v[0], v[1], v[2])
        model_part.CreateNewNode(4, v_2[0], v_2[1], v_2[2])
        model_part.CreateNewNode(5, g[0], g[1], g[2])
        elem_1 = model_part.CreateNewElement('Quad4NodeElement', 1, [1, 2, 3, 4], [5])
        true_area = elem_1.TrueArea()
        area_gauss_1 = elem_1.Area()
        # print("true_area: ", true_area)
        # print("area: ", area)
        ref_area = np.linalg.norm(np.cross(v_1,v_2))
        self.assertEqualTolerance(true_area, ref_area, 1e-10)
        self.assertEqualTolerance(area_gauss_1, ref_area*0.5, 1e-10)

    def test_element_average_calculation_with_vertices_as_nodes(self):
        g1 = 0.0
        g2 = (1 / math.sqrt(3))
        model = data_structure.Model()
        model_part = model.CreateModelPart('test_model_part')
        pressure_var = vars(data_structure)['PRESSURE']
        traction_var = vars(data_structure)['TRACTION']
        model_part.AddNodalSolutionStepVariable(pressure_var)
        model_part.AddNodalSolutionStepVariable(traction_var)
        model_part.CreateNewNode(1, -1.0, -1.0, 0.0)
        model_part.CreateNewNode(2, 1.0, -1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(4, -1.0, 1.0, 0.0)
        gauss_node_1 = model_part.CreateNewNode(5, g2, g2, 0.0)
        gauss_node_2 = model_part.CreateNewNode(6, -1 * g2, g2, 0.0)
        gauss_node_3 = model_part.CreateNewNode(7, -1 * g2, -1 * g2, 0.0)
        gauss_node_4 = model_part.CreateNewNode(8, g2, -1 * g2, 0.0)
        gauss_node_0 = model_part.CreateNewNode(9, g1, g1, 0.0)

        ## Applying linear pressure field whose avg is 0.5
        p_1 = 1.0
        p_2 = 0.0
        gauss_node_1.SetSolutionStepValue(pressure_var,0, p_1)
        gauss_node_2.SetSolutionStepValue(pressure_var,0, p_2)
        gauss_node_3.SetSolutionStepValue(pressure_var,0, p_2)
        gauss_node_4.SetSolutionStepValue(pressure_var,0, p_1)
        gauss_node_0.SetSolutionStepValue(pressure_var, 0,0.5)


        ## Applying linear traction field whose avg is [0.5,0.25,0.0]
        t_1 = [1.0, 0.5, 1.0]
        t_2 = [0.0, 0.0, -1.0]
        gauss_node_1.SetSolutionStepValue(traction_var, 0, t_1)
        gauss_node_2.SetSolutionStepValue(traction_var, 0, t_2)
        gauss_node_3.SetSolutionStepValue(traction_var, 0,  t_2)
        gauss_node_4.SetSolutionStepValue(traction_var, 0, t_1)
        gauss_node_0.SetSolutionStepValue(traction_var, 0, [0.5,0.25,0.0])
        elem_1 = model_part.CreateNewElement('Quad4NodeElement', 1, [1, 2, 3, 4], [5,6,7,8])
        elem_2 = model_part.CreateNewElement('Quad4NodeElement', 2, [1, 2, 3, 4], [9])
        avg_pres_gauss_4 = elem_1.CalculateElementAverage(pressure_var)
        avg_pres_gauss_1 = elem_2.CalculateElementAverage(pressure_var)
        ref_pres_avg =  0.5*(p_1 + p_2)
        # print('ref t_avg: ', ref_pres_avg)
        # print(avg_pres_gauss_4)
        # print(avg_pres_gauss_1)
        self.assertEqualTolerance(avg_pres_gauss_4, ref_pres_avg, 1e-10)
        self.assertEqualTolerance(avg_pres_gauss_1, ref_pres_avg, 1e-10)

        avg_trac_gauss_4 = elem_1.CalculateElementAverage(traction_var)
        avg_trac_gauss_1 = elem_2.CalculateElementAverage(traction_var)
        ref_trac_avg = 0.5*(np.array(t_1)+np.array(t_2))
        # print('ref t_avg: ', ref_trac_avg)
        # print(avg_trac_gauss_4)
        # print(avg_trac_gauss_1)
        for i in range(ref_trac_avg.size):
            self.assertEqualTolerance(avg_trac_gauss_4[i], ref_trac_avg[i], 1e-10)
            self.assertEqualTolerance(avg_trac_gauss_1[i], ref_trac_avg[i], 1e-10)

    def test_element_average_calculation_with_vertices_as_points(self):
        g1 = 0.0
        g2 = (1 / math.sqrt(3))
        model = data_structure.Model()
        model_part = model.CreateModelPart('test_model_part')
        pressure_var = vars(data_structure)['PRESSURE']
        traction_var = vars(data_structure)['TRACTION']
        model_part.AddNodalSolutionStepVariable(pressure_var)
        model_part.AddNodalSolutionStepVariable(traction_var)
        model_part.CreateNewPoint(1, -1.0, -1.0, 0.0)
        model_part.CreateNewPoint(2, 1.0, -1.0, 0.0)
        model_part.CreateNewPoint(3, 1.0, 1.0, 0.0)
        model_part.CreateNewPoint(4, -1.0, 1.0, 0.0)
        gauss_node_1 = model_part.CreateNewNode(5, g2, g2, 0.0)
        gauss_node_2 = model_part.CreateNewNode(6, -1 * g2, g2, 0.0)
        gauss_node_3 = model_part.CreateNewNode(7, -1 * g2, -1 * g2, 0.0)
        gauss_node_4 = model_part.CreateNewNode(8, g2, -1 * g2, 0.0)
        gauss_node_0 = model_part.CreateNewNode(9, g1, g1, 0.0)

        ## Applying linear pressure field whose avg is 0.5
        p_1 = 1.0
        p_2 = 0.0
        gauss_node_1.SetSolutionStepValue(pressure_var,0, p_1)
        gauss_node_2.SetSolutionStepValue(pressure_var,0, p_2)
        gauss_node_3.SetSolutionStepValue(pressure_var,0, p_2)
        gauss_node_4.SetSolutionStepValue(pressure_var,0, p_1)
        gauss_node_0.SetSolutionStepValue(pressure_var, 0,0.5)


        ## Applying linear traction field whose avg is [0.5,0.25,0.0]
        t_1 = [1.0, 0.5, 1.0]
        t_2 = [0.0, 0.0, -1.0]
        gauss_node_1.SetSolutionStepValue(traction_var, 0, t_1)
        gauss_node_2.SetSolutionStepValue(traction_var, 0, t_2)
        gauss_node_3.SetSolutionStepValue(traction_var, 0,  t_2)
        gauss_node_4.SetSolutionStepValue(traction_var, 0, t_1)
        gauss_node_0.SetSolutionStepValue(traction_var, 0, [0.5,0.25,0.0])
        elem_1 = model_part.CreateNewElement('Quad4NodeElement', 1, [1, 2, 3, 4], [5,6,7,8])
        elem_2 = model_part.CreateNewElement('Quad4NodeElement', 2, [1, 2, 3, 4], [9])
        avg_pres_gauss_4 = elem_1.CalculateElementAverage(pressure_var)
        avg_pres_gauss_1 = elem_2.CalculateElementAverage(pressure_var)
        ref_pres_avg =  0.5*(p_1 + p_2)
        # print('ref t_avg: ', ref_pres_avg)
        # print(avg_pres_gauss_4)
        # print(avg_pres_gauss_1)
        self.assertEqualTolerance(avg_pres_gauss_4, ref_pres_avg, 1e-10)
        self.assertEqualTolerance(avg_pres_gauss_1, ref_pres_avg, 1e-10)

        avg_trac_gauss_4 = elem_1.CalculateElementAverage(traction_var)
        avg_trac_gauss_1 = elem_2.CalculateElementAverage(traction_var)
        ref_trac_avg = 0.5*(np.array(t_1)+np.array(t_2))
        # print('ref t_avg: ', ref_trac_avg)
        # print(avg_trac_gauss_4)
        # print(avg_trac_gauss_1)
        for i in range(ref_trac_avg.size):
            self.assertEqualTolerance(avg_trac_gauss_4[i], ref_trac_avg[i], 1e-10)
            self.assertEqualTolerance(avg_trac_gauss_1[i], ref_trac_avg[i], 1e-10)








if __name__ == '__main__':
        KratosUnittest.main()
