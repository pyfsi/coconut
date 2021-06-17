import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import os
import time
import pandas as pd
import numpy as np


class StructuralMechanicsWrapper(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        self.interfaces = [elem.GetString() for elem in project_parameters["interface_sub_model_parts_list"]]
        super(StructuralMechanicsAnalysis, self).__init__(model, project_parameters)
        self.coupling_iteration = None
        self.initial_displacement =  {}

    def Initialize(self):
        super(StructuralMechanicsAnalysis, self).Initialize()

        for sub_mp_name in self.interfaces:
            file_name = f'{sub_mp_name}_nodes.csv'
            with open(file_name, 'w') as f:
                f.write('node_id, x0, y0, z0\n')
            sub_model_part = self.GetSubModelPart(sub_mp_name)
            init_disp = self.GetInitialNodalDisplacement(sub_model_part_name=sub_mp_name)
            with open(file_name, 'a') as f:
                for i, node in enumerate(sub_model_part.Nodes):
                    f.write(f'{node.Id}, {node.X0 + init_disp[i, 0]}, {node.Y0 + init_disp[i, 1]}, {node.Z0 + init_disp[i, 2]}\n')

    def InitializeSolutionStep(self):
        self.time = self._GetSolver().AdvanceInTime(self.time)
        super(StructuralMechanicsAnalysis, self).InitializeSolutionStep()
        self._GetSolver().Predict()
        self.coupling_iteration = 0

    def SolveSolutionStep(self):
        self.coupling_iteration += 1
        KratosMultiphysics.Logger.Print(f'Coupling iteration: {self.coupling_iteration}')
        self.InputData()
        self._GetSolver().SolveSolutionStep()
        KratosMultiphysics.Logger.Print(f'Coupling iteration {self.coupling_iteration} end')
        self.OutputData()


    def OutputData(self):

        for sub_model_part_name in self.interfaces:
            init_disp = self.GetInitialNodalDisplacement(sub_model_part_name=sub_model_part_name)
            full_sub_model_part_name = "Structure." + sub_model_part_name
            if self.model["Structure"].HasSubModelPart(sub_model_part_name):
                sub_model_part = self.model[full_sub_model_part_name]
                file_name = f'{sub_model_part_name}_displacement.csv'
                with open(file_name, 'w') as f:
                    f.write('node_id, displacement_x, displacement_y, displacement_z\n')

                for i, node in enumerate(sub_model_part.Nodes):
                    with open(file_name, 'a') as f:
                        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                        f.write(str(node.Id) + ', ' + str(disp[0]-init_disp[i, 0]) + ', ' + str(disp[1]-init_disp[i, 1]) + ', ' + str(disp[2]-init_disp[i, 2]) + '\n')
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
                        sub_model_part.GetNode(node_id).SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE, 0,
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

    def GetInitialNodalDisplacement(self, sub_model_part_name):
        sub_model_part = self.GetSubModelPart(sub_model_part_name)
        file_name = f'{sub_model_part_name}_init_displacement.csv'
        try:
            init_disp = np.loadtxt(file_name, comments='#')
        except IOError:
            init_disp = np.zeros((sub_model_part.NumberOfNodes(), 3))

        return init_disp

    def GetSubModelPart(self, sub_model_part_name):
        full_sub_model_part_name = "Structure." + sub_model_part_name
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
    model = KratosMultiphysics.Model()
    with open(parameter_file_name, 'r') as parameter_file:
        kratos_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    str_wrapper = StructuralMechanicsWrapper(model, kratos_parameters)
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
