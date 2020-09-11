from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface
from coconut.coupling_components import tools
import numpy as np

from coconut.coupling_components import output_interface_vtk_process

def TestOutputInterfaceVtkProcess():

    parameter_file_name = 'test_interface_output.json'

    with open(parameter_file_name, 'r') as parameter_file:
        parameters = data_structure.Parameters(parameter_file.read())

    input_interface_param = parameters["interface_input"]
    output_interface_param = parameters["interface_output"]

    variable_pres = vars(data_structure)["PRESSURE"]
    variable_trac = vars(data_structure)["TRACTION"]
    # variable_disp = vars(data_structure)["DISPLACEMENT"]

    input_interface_names = input_interface_param.keys()
    #output_interface_names = output_interface_param.keys()



    model = data_structure.Model()
    model_part = model.CreateModelPart(input_interface_param.keys()[0])

    ## add variables to model_part

    model_part.AddNodalSolutionStepVariable(variable_pres)
    model_part.AddNodalSolutionStepVariable(variable_trac)

    # create nodes
    model_part.CreateNewPoint(1, 0.00000, 1.00000, 0.00000)
    model_part.CreateNewPoint(2, 0.00000, 0.50000, 0.00000)
    model_part.CreateNewPoint(3, 0.50000, 1.00000, 0.00000)
    model_part.CreateNewPoint(4, 0.50000, 0.50000, 0.00000)
    model_part.CreateNewPoint(5, 0.00000, 0.00000, 0.00000)
    model_part.CreateNewPoint(6, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewPoint(7, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewPoint(8, 0.50000, 0.00000, 0.00000)
    model_part.CreateNewPoint(9, 1.00000, 0.00000, 0.00000)
    model_part.CreateNewPoint(10, 1.50000, 1.00000, 0.00000)
    model_part.CreateNewPoint(11, 1.50000, 0.50000, 0.00000)
    model_part.CreateNewPoint(12, 1.50000, 0.00000, 0.00000)
    model_part.CreateNewPoint(13, 2.00000, 1.00000, 0.00000)
    model_part.CreateNewPoint(14, 2.00000, 0.50000, 0.00000)
    model_part.CreateNewPoint(15, 2.00000, 0.00000, 0.00000)
    model_part.CreateNewPoint(16, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewPoint(17, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewPoint(18, 1.00000, 0.00000, 0.00000)

    # create Element
    def CreateGaussNode(point_ids, id):
        gauss_coordinates = np.zeros(3)
        for point_id in point_ids:
            gauss_coordinates += model_part.GetPoint(point_id).Coordinates()
        gauss_coordinates /= 4
        model_part.CreateNewNode(id, gauss_coordinates[0], gauss_coordinates[1], gauss_coordinates[2] )


    element_point_ids = [14, 11, 12, 15]
    gauss_node_id = 19
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 1,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [13, 10, 11, 14]
    gauss_node_id = 20
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 2,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [11, 17, 18, 12]
    gauss_node_id = 21
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 3,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [10, 16, 17, 11]
    gauss_node_id = 22
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 4,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [2, 4, 3, 1]
    gauss_node_id = 23
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 5,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [5, 8, 4, 2]
    gauss_node_id = 24
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 6,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [4, 7, 6, 3]
    gauss_node_id = 25
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 7,
                                element_point_ids, [gauss_node_id])

    element_point_ids = [8, 9, 7, 4]
    gauss_node_id = 26
    CreateGaussNode(element_point_ids, gauss_node_id)
    model_part.CreateNewElement("Quad4NodeElement", 8,
                                element_point_ids, [gauss_node_id])


    for node in model_part.Nodes:
        node.SetSolutionStepValue(variable_pres,0,1000)
        node.SetSolutionStepValue(variable_trac,0, [100, 50, 10])

    interface_input = Interface(model, input_interface_param)

    out_obj = output_interface_vtk_process.OutputInterfaceVtkProcess(interface_input)
    time_step = 0
    iteration = 0
    out_obj.Write(time_step, iteration)

if __name__ == '__main__':
    TestOutputInterfaceVtkProcess()