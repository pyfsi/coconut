import KratosMultiphysics as KM
from KratosMultiphysics import StructuralMechanicsApplication as SM
from structural_mechanics_analysis import StructuralMechanicsAnalysis

import os
import time
import json
import argparse
import pandas as pd


class StructuralMechanicsWrapper:

    def __init__(self, parameter_file_name):

        self.model = KM.Model()
        with open(parameter_file_name, 'r') as parameter_file:
            self.kratos_parameters = KM.Parameters(parameter_file.read())

        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)

        self.interfaces = parameters["interface_sub_model_parts_list"]

        self.structural_analysis = StructuralMechanicsAnalysis(self.model, self.kratos_parameters)

    def Initialize(self):
        self.structural_analysis.Initialize()

        for sub_model_part_name in self.interfaces:

            sub_model_part = self.model[f'Structure.{sub_model_part_name}']
            file_name = f'{sub_model_part_name}_nodes.csv'
            with open(file_name, 'w') as f:
                f.write('node_id, x0, y0, z0\n')

            for node in sub_model_part.Nodes:
                with open(file_name, 'a') as f:
                    f.write(f'{node.Id}, {node.X0}, {node.Y0}, {node.Z0}\n')

    def InitializeSolutionStep(self):
        self.structural_analysis.time = self.structural_analysis.solver.AdvanceInTime(self.structural_analysis.time)
        self.structural_analysis.InitializeSolutionStep()
        self.structural_analysis.solver.Predict()

    def SolveSolutionStep(self):
        self.InputData()
        self.structural_analysis.solver.SolveSolutionStep()
        self.OutputData()

    def FinalizeSolutionStep(self):
        self.structural_analysis.FinalizeSolutionStep()

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
                with open(file_name, 'w') as f:
                    f.write('node_id, displacement_x, displacement_y, displacement_z\n')

                for node in sub_model_part.Nodes:
                    with open(file_name, 'a') as f:
                        disp = node.GetSolutionStepValue(KM.DISPLACEMENT)
                        f.write(str(node.Id) + ', ' + str(disp[0]) + ', ' + str(disp[1]) + ', ' + str(disp[2]) + '\n')


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
                    node_ids = pressure_data.node_id
                    pressure = pressure_data.pressure
                    for i, node_id in enumerate(node_ids):
                        sub_model_part.GetNode(node_id).SetSolutionStepValue(KM.POSITIVE_FACE_PRESSURE, 0,
                                                                             -1 * pressure[i])

                file_name_sl = f'{sub_model_part_name}_surface_load.csv'
                if os.path.isfile(file_name_sl):
                    surface_load_data = pd.read_csv(file_name_sl, skipinitialspace=True)
                    node_ids = surface_load_data.node_id
                    surface_load_x = surface_load_data.surface_load_x
                    surface_load_y = surface_load_data.surface_load_y
                    surface_load_z = surface_load_data.surface_load_z
                    for i, node_id in enumerate(node_ids):
                        sub_model_part.GetNode(node_id).SetSolutionStepValue(SM.SURFACE_LOAD, 0,
                                                                             [surface_load_x[i], surface_load_y[i],
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

    while (True):
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
