from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface
from coconut.coupling_components import tools

import os
from os.path import join
import json
import time
import pandas as pd
from subprocess import Popen


def Create(parameters):
    return SolverWrapperKratosStructure6_0(parameters)


class SolverWrapperKratosStructure6_0(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.working_directory = join(os.getcwd(), self.settings["working_directory"].GetString())
        delta_t = self.settings["delta_t"].GetDouble()
        timestep_start = self.settings["timestep_start"].GetDouble()

        input_file_name = join(self.working_directory, self.settings["input_file"].GetString())

        with open(input_file_name, "r") as parameter_file:
            kratos_parameters = json.load(parameter_file)

        kratos_parameters["problem_data"]["start_time"] = timestep_start
        kratos_parameters["problem_data"]["time_step"] = delta_t

        with open(os.path.join(self.working_directory, input_file_name), 'w') as f:
            json.dump(kratos_parameters, f, indent=4)

        self.CheckInterface()

        self.variable_pres = vars(data_structure)["PRESSURE"]
        self.variable_trac = vars(data_structure)["TRACTION"]
        self.variable_disp = vars(data_structure)["DISPLACEMENT"]
        self.model = data_structure.Model()

        input_interface_names = self.settings["interface_input"].keys()
        output_interface_names = self.settings["interface_output"].keys()

        rel_path = 'coconut/coupling_components/solver_wrappers/kratos/run_kratos_structural.py'
        python_file = os.path.join(os.environ['COCONUT_PATH'], rel_path)
        log_file = os.path.join(self.working_directory, 'log')
        self.kratos_process = Popen(f'python3 {python_file} {input_file_name} &> {log_file}', shell=True,
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

        self.iteration = 0

        self.send_message('next')
        self.wait_message('next_ready')

    def SolveSolutionStep(self, interface_input):

        self.iteration += 1

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
        for name in input_interface_list:
            if '_input' not in name:
                raise RuntimeError(
                    '<sub_mp_name> in the list of "sub_model_parts_from_str_solver" in json file should have entry as <sub_mp_name>_input')

        for name in output_interface_list:
            if '_output' not in name:
                raise RuntimeError(
                    '<sub_mp_name> in the list of "sub_model_parts_from_str_solver" in json file should have entry as <sub_mp_name>_output')

        sub_model_part_list = self.settings["sub_model_parts_from_str_solver"].list()
        sub_model_part_list = [elem.GetString() for elem in sub_model_part_list]
        sub_model_part_list_interface = [name.replace('_input', '') for name in input_interface_list]
        sub_model_part_list_interface += [name.replace('_output', '') for name in output_interface_list]
        if not set(sub_model_part_list) == set(sub_model_part_list_interface):
            raise RuntimeError(
                '<name> in the list of "sub_model_parts_from_str_solver" in json file should have entry as <name_input> in "interface_input" or <name_output> in "interface_output" ')

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
