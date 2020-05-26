from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface
from coconut.coupling_components import tools

from coconut.coupling_components import output_vtk_process

def TestOutputVtk():

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
    model_part.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
    model_part.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
    model_part.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
    model_part.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
    model_part.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
    model_part.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
    model_part.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
    model_part.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
    model_part.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
    model_part.CreateNewNode(18, 1.00000, 0.00000, 0.00000)

    # create Element
    model_part.CreateNewElement("Quad", 1,
                                [14, 11, 12, 15], 0)
    model_part.CreateNewElement("Quad", 2,
                                [13, 10, 11, 14], 0)
    model_part.CreateNewElement("Quad", 3,
                                [11, 17, 18, 12],0)
    model_part.CreateNewElement("Quad", 4,
                                [10, 16, 17, 11], 0)
    model_part.CreateNewElement("Quad", 5,
                                [2, 4, 3, 1], 0)
    model_part.CreateNewElement("Quad", 6, [5, 8, 4, 2],0)
    model_part.CreateNewElement("Quad", 7, [4, 7, 6, 3],0)
    model_part.CreateNewElement("Quad", 8, [8, 9, 7, 4],0)


    for node in model_part.Nodes:
        node.SetSolutionStepValue(variable_pres,0,1000)
        node.SetSolutionStepValue(variable_trac,0, [100, 50, 10])

    interface_input = Interface(model, input_interface_param)

    out_obj = output_vtk_process.OutputInterfaceProcess(interface_input)
    out_obj.Write()

if __name__ == '__main__':
    TestOutputVtk()





