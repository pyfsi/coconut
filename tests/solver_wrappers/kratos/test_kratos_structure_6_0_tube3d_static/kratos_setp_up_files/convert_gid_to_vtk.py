import numpy as np
import glob
import json

geom_type_dict = {4: 9}


def GetLinesBetweenTwoStrings(lines, start_string, end_string, index):
    output_lines = []
    add_lines = False

    for line in lines[index:]:

        if end_string in line and add_lines:
            add_lines = False
            index += 1
            break

        if (add_lines):
            output_lines.append(line)

        if start_string in line:
            add_lines = True

        index += 1

    return output_lines, index;


def GetIndexOfString(lines, search_string):
    for index, line in enumerate(lines):
        if search_string in line:
            return index
    print(f'{search_string} not found!')
    return -1


class GidToVtkConverter:

    def __init__(self, mesh_file_name, result_file_name, output_vtk_file_name):

        with open(mesh_file_name, 'r') as mesh_file:
            self.mesh_lines = mesh_file.readlines()

        with open(result_file_name, 'r') as result_file:
            self.result_lines = result_file.readlines()

        found_quad = False

        for line in self.mesh_lines:
            if 'Quadrilateral3D4' in line:
                found_quad = True

        if not found_quad:
            raise Exception('Other elements not implmented yet')

        self.coordinates = []
        self.element_conns = []
        self.displacement_results = []
        self.pressure_results = []
        self.traction_results = []
        self.variable_values_dict = {'PRESSURE': self.pressure_results,
                                     'TRACTION': self.traction_results,
                                     'DISPLACEMENT': self.displacement_results}

        self.output_vtk_file_name = output_vtk_file_name

        with open(self.output_vtk_file_name, 'w') as output_vtk_file:
            output_vtk_file.write('# vtk DataFile Version 4.0\n')
            output_vtk_file.write('vtk output\n')
            output_vtk_file.write('ASCII\n')
            output_vtk_file.write('DATASET UNSTRUCTURED_GRID\n')

    def ReadNodes(self):
        node_lines, _ = GetLinesBetweenTwoStrings(self.mesh_lines, 'Coordinates', 'End Coordinates', 0)

        for line in node_lines:
            node_coord_list = [float(i) for i in line.split(' ')]
            self.coordinates.append(node_coord_list[1:])

    def ReadElements(self):
        elem_lines, _ = GetLinesBetweenTwoStrings(self.mesh_lines, 'Elements', 'End Elements', 0)

        for line in elem_lines:
            conn_list = [int(i) - 1 for i in line.split(' ')]
            self.element_conns.append(conn_list[1:-1])

    def ReadResults(self):

        index = GetIndexOfString(self.result_lines, 'DISPLACEMENT')

        if not index == -1:
            result_lines, _ = GetLinesBetweenTwoStrings(self.result_lines, 'Values', 'End Values', index)
            for line in result_lines:
                self.displacement_results.append(
                    [float(line.split(' ')[1]), float(line.split(' ')[2]), float(line.split(' ')[3])])

        index = GetIndexOfString(self.result_lines, 'SURFACE_LOAD')
        if not index == -1:
            result_lines, _ = GetLinesBetweenTwoStrings(self.result_lines, 'Values', 'End Values', index)
            for line in result_lines:
                self.traction_results.append(
                    [float(line.split(' ')[1]), float(line.split(' ')[2]), float(line.split(' ')[3])])
        index = GetIndexOfString(self.result_lines, 'POSITIVE_FACE_PRESSURE')
        if not index == -1:
            result_lines, _ = GetLinesBetweenTwoStrings(self.result_lines, 'Values', 'End Values', index)
            for line in result_lines:
                self.pressure_results.append([float(line.split(' ')[1])])

    def WriteNodes(self):
        with open(self.output_vtk_file_name, 'a') as output_vtk_file:
            output_vtk_file.write(f'POINTS {len(self.coordinates)} float\n')

            for node in self.coordinates:
                output_vtk_file.write(f'{node[0]} {node[1]} {node[2]}\n')

    def DetermineVtkCellListSize(self):
        vtk_cell_list_size = 0
        for elem_conn in self.element_conns:
            vtk_cell_list_size += len(elem_conn) + 1
        return vtk_cell_list_size

    def WriteElements(self):
        cell_list_size = self.DetermineVtkCellListSize()
        with open(self.output_vtk_file_name, 'a') as output_vtk_file:
            output_vtk_file.write(f'\nCELLS {len(self.element_conns)}  {cell_list_size}\n')

            for elem_conn in self.element_conns:
                output_vtk_file.write(f'{len(elem_conn)} {elem_conn[0]} {elem_conn[1]} {elem_conn[2]} {elem_conn[3]}\n')

            output_vtk_file.write(f'\nCELL_TYPES {len(self.element_conns)}\n')

            for elem_conn in self.element_conns:
                output_vtk_file.write(f'{geom_type_dict[len(elem_conn)]}\n')

    def WriteNodalResults(self):

        node_variables = []
        if self.pressure_results:
            node_variables.append({'variable_name': 'PRESSURE', 'nr_components': 1})
        if self.traction_results:
            node_variables.append({'variable_name': 'TRACTION', 'nr_components': 3})
        if self.displacement_results:
            node_variables.append({'variable_name': 'DISPLACEMENT', 'nr_components': 3})

        with open(self.output_vtk_file_name, 'a') as output_vtk_file:
            output_vtk_file.write(f'\nPOINT_DATA {len(self.coordinates)}\n')
            output_vtk_file.write(f'FIELD FieldData {len(node_variables)}\n')

            for node_variable in node_variables:
                output_vtk_file.write(f'{node_variable["variable_name"]} {node_variable["nr_components"]} {len(
                    self.coordinates)}  float\n')
                values = self.variable_values_dict[node_variable['variable_name']]
                for value in values:
                    for i in range(0, node_variable['nr_components']):
                        output_vtk_file.write(f'{value[i]} ')

                    output_vtk_file.write('\n')

    def Execute(self):
        self.ReadNodes()
        self.ReadElements()
        self.ReadResults()
        self.WriteNodes()
        self.WriteElements()
        self.WriteNodalResults()


if __name__ == '__main__':

    with open('ProjectParameters.json') as param_file:
        param = json.load(param_file)

    problem_name = param["problem_data"]["problem_name"]
    for file_path in glob.iglob('*.res'):
        step = int(file_path.split('.')[0].split('_')[-1])
        msh_file = f'{problem_name}_{step}.post.msh'
        res_file = f'{problem_name}_{step}.post.res'
        vtk_file = f'{problem_name}_{step}.vtk'
        gid_converter = GidToVtkConverter(msh_file, res_file, vtk_file)
        gid_converter.Execute()
