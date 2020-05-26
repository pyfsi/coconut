import numpy as np
import glob
import json
from coconut import data_structure
import os




class OutputInterfaceProcess:

    def __init__(self, interface, path = './'):


        self.interface = interface
        self.model_part_variables_dict = {}
        self.step = 0
        self.path = path

        for model_part_name, variable_names in interface.model_parts_variables:
            variable_list = []
            for variable_name in variable_names.list():
                variable_list.append(vars(data_structure)[variable_name.GetString()])
            self.model_part_variables_dict[model_part_name] = variable_list
        self.geom_type_dict = {4: 9}



    def WriteHeader(self):

        for model_part_name, _  in self.model_part_variables_dict.items():


            output_file_name =f'{model_part_name}_{self.step}.vtk'
            output_file_name = os.path.join(self.path, output_file_name)

        with open(output_file_name, 'w') as output_vtk_file:
            output_vtk_file.write('# vtk DataFile Version 4.0\n')
            output_vtk_file.write('vtk output\n')
            output_vtk_file.write('ASCII\n')
            output_vtk_file.write('DATASET UNSTRUCTURED_GRID\n')



    def Write(self):

        self.WriteHeader()
        self.WriteNodes()
        self.WriteElements()
        self.WriteNodalResults()
        #self.ElementResults()
        self.step += 1



    def WriteNodes(self):
        for model_part_name, _ in self.model_part_variables_dict.items():

            output_file_name = f'{model_part_name}_{self.step}.vtk'
            output_file_name = os.path.join(self.path, output_file_name)
            model_part = self.interface.model.GetModelPart(model_part_name)

            with open(output_file_name, 'a') as output_vtk_file:
                output_vtk_file.write(f'POINTS {model_part.NumberOfNodes()} float\n')

                for node in model_part.Nodes:
                    output_vtk_file.write(f'{node.X0} {node.Y0} {node.Z0}\n')


    def DetermineVtkCellListSize(self, model_part):
        vtk_cell_list_size = 0
        for elem in model_part.Elements:
            vtk_cell_list_size += elem.NumberOfNodes()+1
        return vtk_cell_list_size




    def WriteElements(self):

        for model_part_name, _ in self.model_part_variables_dict.items():

            output_file_name = f'{model_part_name}_{self.step}.vtk'
            output_file_name = os.path.join(self.path, output_file_name)
            model_part = self.interface.model.GetModelPart(model_part_name)
            cell_list_size = self.DetermineVtkCellListSize(model_part)

            elements = model_part.Elements

            with open(output_file_name, 'a') as output_vtk_file:
                output_vtk_file.write(f'\nCELLS {model_part.NumberOfElements()}  {cell_list_size}\n')

                for elem in model_part.Elements:
                    output_vtk_file.write(f'{elem.NumberOfNodes()} ')

                    nodes = elem.GetNodes()
                    for node in nodes:
                        output_vtk_file.write(f'{node.Id-1} ')
                    output_vtk_file.write('\n')

                output_vtk_file.write(f'\nCELL_TYPES {model_part.NumberOfElements()}\n')

                for elem in model_part.Elements:
                    output_vtk_file.write(f'{self.geom_type_dict[elem.NumberOfNodes()]}\n')



    def WriteNodalResults(self):

        for model_part_name, variable_list in self.model_part_variables_dict.items():
            output_file_name = f'{model_part_name}_{self.step}.vtk'
            output_file_name = os.path.join(self.path, output_file_name)
            model_part = self.interface.model.GetModelPart(model_part_name)
            cell_list_size = self.DetermineVtkCellListSize(model_part)
            variable_components_dict = {"Array" : 3, "Double" : 1}

            with open(output_file_name, 'a') as output_vtk_file:
                    output_vtk_file.write(f'\nPOINT_DATA {model_part.NumberOfNodes()}\n')
                    output_vtk_file.write(f'FIELD FieldData {len(variable_list)}\n')

                    for variable in variable_list:
                        output_vtk_file.write(f'{variable.Name()} {variable_components_dict[variable.Type()]} {model_part.NumberOfNodes()}  float\n')

                        for node in model_part.Nodes:
                                value = node.GetSolutionStepValue(variable)
                                if variable.Type() == "Double":
                                    output_vtk_file.write(f'{value} ')
                                elif variable.Type() ==  "Array":
                                    for component in value:
                                        output_vtk_file.write(f'{component} ')

                                output_vtk_file.write('\n')

    def ElementResults(self):
        pass






