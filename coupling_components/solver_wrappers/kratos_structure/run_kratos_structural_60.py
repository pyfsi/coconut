import KratosMultiphysics as KM
from KratosMultiphysics import StructuralMechanicsApplication as SM
from structural_mechanics_analysis import StructuralMechanicsAnalysis

import os
import time
import json
import pandas as pd
import numpy as np


class StructuralMechanicsWrapper:

    def __init__(self, parameter_file_name):

        self.model = KM.Model()
        with open(parameter_file_name, 'r') as parameter_file:
            self.kratos_parameters = KM.Parameters(parameter_file.read())

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        self.coupling_iteration = None
        self.interfaces = parameters["interface_sub_model_parts_list"]

        self.structural_analysis = StructuralMechanicsAnalysis(self.model, self.kratos_parameters)

    def Initialize(self):
        self.structural_analysis.Initialize()

        for sub_model_part_name in self.interfaces:
            sub_model_part = self.model[f'Structure.{sub_model_part_name}']
            file_name_nodes = f'{sub_model_part_name}_nodes.csv'
            file_name_cond = f'{sub_model_part_name}_conditions.csv'
            node_ids = np.array([node.Id for node in sub_model_part.Nodes])
            cond_ids = np.array([cond.Id for cond in sub_model_part.Conditions])
            cond_centres = np.array(
                [[cond.GetGeometry().Center().X, cond.GetGeometry().Center().Y, cond.GetGeometry().Center().Z] for cond
                 in sub_model_part.Conditions])
            node_coords0 = np.array([[node.X0, node.Y0, node.Z0] for node in sub_model_part.Nodes])
            node_coords0_df = pd.DataFrame(
                {'node_id': node_ids, 'x0': node_coords0[:, 0], 'y0': node_coords0[:, 1], 'z0': node_coords0[:, 2]})
            node_coords0_df.to_csv(file_name_nodes, index=False)
            cond_coords_df = pd.DataFrame(
                {'cond_id': cond_ids, 'centre_x0': cond_centres[:, 0], 'centre_y0': cond_centres[:, 1],
                 'centre_z0': cond_centres[:, 2]})
            cond_coords_df.to_csv(file_name_cond, index=False)

    def InitializeSolutionStep(self):
        self.structural_analysis.time = self.structural_analysis.solver.AdvanceInTime(self.structural_analysis.time)
        self.structural_analysis.InitializeSolutionStep()
        self.structural_analysis.solver.Predict()
        self.coupling_iteration = 0

    def SolveSolutionStep(self):
        self.coupling_iteration += 1
        KM.Logger.Print(f'Coupling iteration: {self.coupling_iteration}')
        self.InputData()
        self.structural_analysis.solver.SolveSolutionStep()
        KM.Logger.Print(f'Coupling iteration {self.coupling_iteration} end')
        self.OutputData()

    def FinalizeSolutionStep(self):
        self.structural_analysis.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def OutputSolutionStep(self):
        self.structural_analysis.OutputSolutionStep()

    def Finalize(self):
        self.structural_analysis.Finalize()

    def OutputData(self):

        for sub_model_part_name in self.interfaces:
            full_sub_model_part_name = "Structure." + sub_model_part_name
            if self.model["Structure"].HasSubModelPart(sub_model_part_name):
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
                raise Exception(f"{sub_model_part_name} not present in the Kratos model.")

    def InputData(self):
        for sub_model_part_name in self.interfaces:
            full_sub_model_part_name = "Structure." + sub_model_part_name
            if self.model["Structure"].HasSubModelPart(sub_model_part_name):
                sub_model_part = self.model[full_sub_model_part_name]
                file_name_pr = f'{sub_model_part_name}_pressure.csv'
                if os.path.isfile(file_name_pr):
                    pressure_data = pd.read_csv(file_name_pr, skipinitialspace=True)
                    cond_ids = pressure_data.cond_id
                    pressure = pressure_data.pressure
                    for i, cond_id in enumerate(cond_ids):
                        sub_model_part.GetCondition(cond_id).SetValue(KM.POSITIVE_FACE_PRESSURE, -1 * pressure[i])

                file_name_sl = f'{sub_model_part_name}_surface_load.csv'
                if os.path.isfile(file_name_sl):
                    surface_load_data = pd.read_csv(file_name_sl, skipinitialspace=True)
                    cond_ids = surface_load_data.cond_id
                    surface_load_x = surface_load_data.surface_load_x
                    surface_load_y = surface_load_data.surface_load_y
                    surface_load_z = surface_load_data.surface_load_z
                    for i, cond_id in enumerate(cond_ids):
                        sub_model_part.GetNode(cond_id).SetValue(SM.SURFACE_LOAD, [surface_load_x[i], surface_load_y[i],
                                                                                   surface_load_z[i]])

            else:
                raise Exception(f"{sub_model_part_name} not present in the Kratos model.")


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

    str_wrapper = StructuralMechanicsWrapper(parameter_file_name)
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
            os.remove('save.coco')
            open(os.path.join('save_ready.coco'), 'w').close()

        if os.path.isfile('stop.coco'):
            str_wrapper.Finalize()
            os.remove('stop.coco')
            open('stop_ready.coco', 'w').close()
            break
