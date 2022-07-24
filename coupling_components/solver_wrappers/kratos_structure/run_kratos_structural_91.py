import KratosMultiphysics as KM
from KratosMultiphysics import StructuralMechanicsApplication as SM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import os
import time
import pandas as pd
import numpy as np


class StructuralMechanicsWrapper91(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        self.interfaces = [elem.GetString() for elem in project_parameters['interface_sub_model_parts_list']]
        super().__init__(model, project_parameters)
        self.dimensions = project_parameters['solver_settings']['domain_size'].GetInt()
        self.coupling_iteration = None
        self.pressure_directions = [direc.GetInt() for direc in project_parameters['pressure_directions']]

    def Initialize(self):
        super().Initialize()

        for sub_mp_name in self.interfaces:
            file_name_nodes = f'{sub_mp_name}_nodes.csv'
            sub_model_part = self.GetSubModelPart(sub_mp_name)
            node_ids = np.array([node.Id for node in sub_model_part.Nodes])
            node_coords = np.array([[node.X0, node.Y0, node.Z0] for node in sub_model_part.Nodes])
            node_coords_df = pd.DataFrame(
                {'node_id': node_ids, 'x0': node_coords[:, 0], 'y0': node_coords[:, 1], 'z0': node_coords[:, 2]})
            node_coords_df.to_csv(file_name_nodes, index=False)

    def InitializeSolutionStep(self):
        self.time = self._GetSolver().AdvanceInTime(self.time)
        super().InitializeSolutionStep()
        self._GetSolver().Predict()
        self.coupling_iteration = 0

    def SolveSolutionStep(self):
        self.coupling_iteration += 1
        KM.Logger.Print(f'Coupling iteration: {self.coupling_iteration}')
        self.InputData()
        self._GetSolver().SolveSolutionStep()
        KM.Logger.Print(f'Coupling iteration {self.coupling_iteration} end')
        self.OutputData()

    def OutputData(self):
        for sub_model_part_name in self.interfaces:
            full_sub_model_part_name = 'Structure.' + sub_model_part_name
            if self.model['Structure'].HasSubModelPart(sub_model_part_name):
                sub_model_part = self.model[full_sub_model_part_name]
                file_name = f'{sub_model_part_name}_displacement.csv'
                node_ids = np.array([node.Id for node in sub_model_part.Nodes])
                displacement = np.array(
                    [list(node.GetSolutionStepValue(KM.DISPLACEMENT)) for node in sub_model_part.Nodes])
                disp_df = pd.DataFrame(
                    {'node_id': node_ids, 'displacement_x': displacement[:, 0], 'displacement_y': displacement[:, 1],
                     'displacement_z': displacement[:, 2]})
                disp_df.to_csv(file_name, index=False)
            else:
                raise Exception(f'{sub_model_part_name} not present in the Kratos model.')

    def InputData(self):
        for sub_model_part_name, pressure_direction in zip(self.interfaces, self.pressure_directions):
            full_sub_model_part_name = 'Structure.' + sub_model_part_name
            if self.model['Structure'].HasSubModelPart(sub_model_part_name):
                sub_model_part = self.model[full_sub_model_part_name]
                file_name_pr = f'{sub_model_part_name}_pressure.csv'
                if os.path.isfile(file_name_pr):
                    pressure_data = pd.read_csv(file_name_pr, skipinitialspace=True)
                    node_ids = pressure_data.node_id
                    pressure = pressure_data.pressure
                    for i, node_id in enumerate(node_ids):
                        sub_model_part.GetNode(node_id).SetSolutionStepValue(KM.NEGATIVE_FACE_PRESSURE, 0,
                                                                             pressure_direction * pressure[i])

                file_name_tr = f'{sub_model_part_name}_traction.csv'
                if os.path.isfile(file_name_tr):
                    traction_data = pd.read_csv(file_name_tr, skipinitialspace=True)
                    node_ids = traction_data.node_id
                    traction_x = traction_data.traction_x
                    traction_y = traction_data.traction_y
                    traction_z = traction_data.traction_z
                    kratos_variable = SM.SURFACE_LOAD if self.dimensions == 3 else SM.LINE_LOAD
                    for i, node_id in enumerate(node_ids):
                        sub_model_part.GetNode(node_id).SetSolutionStepValue(kratos_variable, 0,
                                                                             [traction_x[i], traction_y[i],
                                                                              traction_z[i]])
            else:
                raise Exception(f'{sub_model_part_name} not present in the Kratos model.')

    def GetSubModelPart(self, sub_model_part_name):
        full_sub_model_part_name = 'Structure.' + sub_model_part_name
        return self.model[full_sub_model_part_name]


if __name__ == '__main__':
    from sys import argv

    # Check number of command line arguments
    if len(argv) != 2:
        err_msg = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python run_simulation.py <cosim-parameter-file>.json"\n'
        raise Exception(err_msg)

    # Import data structure
    parameter_file_name = argv[1]
    model = KM.Model()
    with open(parameter_file_name, 'r') as parameter_file:
        kratos_parameters = KM.Parameters(parameter_file.read())

    str_wrapper = StructuralMechanicsWrapper91(model, kratos_parameters)
    str_wrapper.Initialize()
    open('start_ready.coco', 'w').close()

    while True:
        time.sleep(0.01)

        if os.path.isfile('next.coco'):
            str_wrapper.InitializeSolutionStep()
            os.remove('next.coco')
            open(os.path.join('next_ready.coco'), 'w').close()

        if os.path.isfile('continue.coco'):
            str_wrapper.SolveSolutionStep()
            os.remove(os.path.join('continue.coco'))
            open(os.path.join('continue_ready.coco'), 'w').close()

        if os.path.isfile('save.coco'):
            str_wrapper.FinalizeSolutionStep()
            str_wrapper.OutputSolutionStep()
            os.remove('save.coco')
            open(os.path.join('save_ready.coco'), 'w').close()

        if os.path.isfile('stop.coco'):
            str_wrapper.Finalize()
            os.remove('stop.coco')
            open('stop_ready.coco', 'w').close()
            break
