import json
import os
import time
from os.path import join
from subprocess import Popen

import pandas as pd
from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface


def Create(parameters):
    return SolverWrapperKratosStructure60(parameters)


class SolverWrapperKratosStructure60(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.working_directory = join(os.getcwd(), self.settings["working_directory"].GetString())
        delta_t = self.settings["delta_t"].GetDouble()
        timestep_start = self.settings["timestep_start"].GetDouble()
        dimensions = self.settings["dimensions"].GetDouble()

        input_file_name = join(self.working_directory, self.settings["input_file"].GetString())

        with open(input_file_name, "r") as parameter_file:
            kratos_parameters = json.load(parameter_file)

        kratos_parameters["problem_data"]["start_time"] = timestep_start
        kratos_parameters["problem_data"]["time_step"] = delta_t
        kratos_parameters["problem_data"]["domain_size"] = dimensions
        interface_sub_model_parts_list = []
        for param in self.settings["kratos_interface_sub_model_parts_list"].list():
            interface_sub_model_parts_list.append(param.GetString())

        kratos_parameters["interface_sub_model_parts_list"] = interface_sub_model_parts_list

        with open(os.path.join(self.working_directory, input_file_name), 'w') as f:
            json.dump(kratos_parameters, f, indent=4)

        self.CheckInterface()

        self.variable_pres = vars(data_structure)["PRESSURE"]
        self.variable_trac = vars(data_structure)["TRACTION"]
        self.variable_disp = vars(data_structure)["DISPLACEMENT"]
        self.model = data_structure.Model()

        input_interface_names = self.settings["interface_input"].keys()
        output_interface_names = self.settings["interface_output"].keys()

        kratos_load_cmd = self.settings["solver_load_cmd"].GetString()
        dir_path = os.path.dirname(os.path.realpath(__file__))
        run_script_file = os.path.join(dir_path, 'run_kratos_structural_60.py')

        self.kratos_process = Popen(f'{kratos_load_cmd} && python3 {run_script_file} {input_file_name} &> log',
                                    shell=True,
                                    cwd=self.working_directory)

        self.wait_message('start_ready')

        for interface_name in input_interface_names:
            model_part = self.model.CreateModelPart(interface_name)
            model_part.AddNodalSolutionStepVariable(self.variable_pres)
            model_part.AddNodalSolutionStepVariable(self.variable_trac)
            model_part_name = interface_name.replace('_input', '')
            file_path = os.path.join(self.working_directory, model_part_name + '_nodes.csv')
            node_data = pd.read_csv(file_path, skipinitialspace=True)
            node_ids = node_data.node_id
            x0 = node_data.x0
            y0 = node_data.y0
            z0 = node_data.z0
            for i, node_id in enumerate(node_ids):
                model_part.CreateNewNode(node_id, x0[i], y0[i], z0[i])

        for interface_name in output_interface_names:
            model_part = self.model.CreateModelPart(interface_name)
            model_part.AddNodalSolutionStepVariable(self.variable_disp)
            model_part_name = interface_name.replace('_output', '')
            file_path = os.path.join(self.working_directory, model_part_name + '_nodes.csv')
            node_data = pd.read_csv(file_path, skipinitialspace=True)
            node_ids = node_data.node_id
            x0 = node_data.x0
            y0 = node_data.y0
            z0 = node_data.z0
            for i, node_id in enumerate(node_ids):
                model_part.CreateNewNode(node_id, x0[i], y0[i], z0[i])

        # # Interfaces
        self.interface_input = Interface(self.model, self.settings["interface_input"])
        self.interface_output = Interface(self.model, self.settings["interface_output"])
        # run time
        self.run_time = 0.0

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.send_message('next')
        self.wait_message('next_ready')

    def SolveSolutionStep(self, interface_input):

        self.interface_input.SetPythonList(interface_input.GetPythonList())
        self.WriteInputData()
        self.send_message('continue')
        self.wait_message('continue_ready')
        self.UpdateInterfaceOutput()
        interface_output = self.GetInterfaceOutput()

        return interface_output

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.send_message('save')
        self.wait_message('save_ready')

    def Finalize(self):
        super().Finalize()
        self.send_message('stop')
        self.wait_message('stop_ready')
        self.remove_all_messages()
        self.kratos_process.kill()

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    def WriteInputData(self):
        input_interface_names = self.settings["interface_input"].keys()

        for interface_name in input_interface_names:
            model_part = self.model.GetModelPart(interface_name)
            sub_model_part_name = interface_name.replace('_input', '')
            file_path_pr = os.path.join(self.working_directory, sub_model_part_name + '_pressure.csv')
            with open(file_path_pr, 'w') as f:
                f.write('node_id, pressure\n')
            file_path_sl = os.path.join(self.working_directory, sub_model_part_name + '_surface_load.csv')
            with open(file_path_sl, 'w') as f:
                f.write('node_id, surface_load_x, surface_load_y, surface_load_z\n')

            for node in model_part.Nodes:
                pressure = node.GetSolutionStepValue(self.variable_pres)
                surface_load = node.GetSolutionStepValue(self.variable_trac)
                with open(file_path_pr, 'a') as f:
                    f.write(str(node.Id) + ', ' + str(pressure) + '\n')

                with open(file_path_sl, 'a') as f:
                    f.write(str(node.Id) + ', ' + str(surface_load[0]) + ', ' + str(surface_load[1]) + ', ' + str(
                        surface_load[2]) + '\n')

    def UpdateInterfaceOutput(self):
        output_interface_names = self.settings["interface_output"].keys()

        for interface_name in output_interface_names:
            model_part = self.model.GetModelPart(interface_name)
            sub_model_part_name = interface_name.replace('_output', '')
            file_path = os.path.join(self.working_directory, sub_model_part_name + '_displacement.csv')
            disp_data = pd.read_csv(file_path, skipinitialspace=True)
            node_ids = disp_data.node_id
            displacement_x = disp_data.displacement_x
            displacement_y = disp_data.displacement_y
            displacement_z = disp_data.displacement_z
            for i, node_id in enumerate(node_ids):
                model_part.GetNode(node_id).SetSolutionStepValue(self.variable_disp, 0,
                                                                 [displacement_x[i], displacement_y[i],
                                                                  displacement_z[i]])

    def CheckInterface(self):
        input_interface_list = self.settings["interface_input"].keys()
        output_interface_list = self.settings["interface_output"].keys()
        sub_model_part_list = self.settings["kratos_interface_sub_model_parts_list"].list()

        for param in sub_model_part_list:
            name = param.GetString()
            if not f'{name}_input' in input_interface_list:
                raise RuntimeError(
                    f'Error in json file: {name}_input not listed in "interface_input": {input_interface_list}.\n. <sub_mp_name> in the "kratos_interface_sub_model_parts_list" in json file should have corresponding <sub_mp_name>_input in "interface_input" list. ')

            if not f'{name}_output' in output_interface_list:
                raise RuntimeError(
                    f'Error in json file: {name}_output not listed in "interface_output": {output_interface_list}.\n. <sub_mp_name> in the "kratos_interface_sub_model_parts_list" in json file should have corresponding <sub_mp_name>_output in "interface_output" list.')

    def send_message(self, message):
        file = join(self.working_directory, message + ".coco")
        open(file, 'w').close()
        return

    def wait_message(self, message):
        file = join(self.working_directory, message + ".coco")
        while not os.path.isfile(file):
            time.sleep(0.01)
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.working_directory, message + ".coco")
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.working_directory):
            if file_name.endswith('.coco'):
                file = join(self.working_directory, file_name)
                os.remove(file)