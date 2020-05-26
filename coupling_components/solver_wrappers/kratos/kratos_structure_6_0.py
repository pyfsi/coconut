from coconut import data_structure
from coconut.coupling_components.component import Component
from coconut.coupling_components.interface import Interface
from coconut.coupling_components import tools

import KratosMultiphysics
from KratosMultiphysics import StructuralMechanicsApplication
import structural_mechanics_analysis

import os
from os.path import join
import numpy as np
import sys
import shutil


def Create(parameters):
    return SolverWrapperKratosStructure6_0(parameters)


class SolverWrapperKratosStructure6_0(Component):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        dir_structure = join(os.getcwd(), self.settings["working_directory"].GetString())
        delta_t = self.settings["delta_t"].GetDouble()
        timestep_start = self.settings["timestep_start"].GetDouble()
        time_integration_method = self.settings["time_integration_method"].GetString()
        time_integration_scheme = self.settings["time_integration_scheme"].GetString()
        input_file_path =  join(dir_structure,self.settings["input_file"].GetString())

        with open(input_file_path, "r") as parameter_file:
            kratos_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        kratos_solver_settings = kratos_parameters["solver_settings"]

        if not delta_t: #### Steady state
            kratos_solver_settings["solver_type"].SetString("static")
            kratos_parameters["problem_data"]["time_step"].SetDouble(1.0)
        else:
            kratos_parameters["problem_data"]["start_time"].SetDouble(timestep_start)
            kratos_solver_settings["solver_type"].SetString("dynamic")
            kratos_parameters["problem_data"]["time_step"].SetDouble(delta_t)

            if kratos_solver_settings.Has("time_integration_method"):
                kratos_solver_settings["time_integration_method"].SetString(time_integration_method)
            else:
                kratos_solver_settings.AddEmptyValue("time_integration_method")
                kratos_solver_settings["time_integration_method"].SetString(time_integration_method)

            if kratos_solver_settings.Has("scheme_type"):
                kratos_solver_settings["scheme_type"].SetString(time_integration_scheme)
            else:
                kratos_solver_settings.AddEmptyValue("scheme_type")
                kratos_solver_settings["scheme_type"].SetString(time_integration_scheme)

        mesh_file_type = kratos_solver_settings["model_import_settings"]["input_type"].GetString()

        if not mesh_file_type == "mdpa":
            raise Exception("Only mesh file of mdpa type is supported")

        mdpa_file_name = kratos_solver_settings["model_import_settings"]["input_filename"].GetString()
        mdpa_file_path = os.path.join(dir_structure, mdpa_file_name)
        kratos_solver_settings["model_import_settings"]["input_filename"].SetString(mdpa_file_path)

        material_file_name = kratos_solver_settings["material_import_settings"]["materials_filename"].GetString()
        material_file_path = os.path.join(dir_structure, material_file_name)
        kratos_solver_settings["material_import_settings"]["materials_filename"].SetString(material_file_path)
        kratos_model = KratosMultiphysics.Model()
        self.structural_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(kratos_model, kratos_parameters)
        self.structural_analysis.Initialize()
        self.kratos_sub_model_part_list = []
        for kratos_sub_mp_name in self.settings["sub_model_parts_from_str_solver"].list():
            input_sub_mp_name = "Structure." + kratos_sub_mp_name.GetString()
            if kratos_model["Structure"].HasSubModelPart(kratos_sub_mp_name.GetString()):
                self.kratos_sub_model_part_list.append(kratos_model[input_sub_mp_name])
            else:
                raise Exception(f"{input_sub_mp_name} not present in the Kratos model.")

        self.CheckInterface()

        self.variable_pres = vars(data_structure)["PRESSURE"]
        self.variable_trac = vars(data_structure)["TRACTION"]
        self.variable_disp = vars(data_structure)["DISPLACEMENT"]
        self.model = data_structure.Model()

        input_interface_names = self.settings["interface_input"].keys()
        output_interface_names = self.settings["interface_output"].keys()

        node_id = 1

        for index, interface_name in enumerate(input_interface_names):
            model_part = self.model.CreateModelPart(interface_name)
            model_part.AddNodalSolutionStepVariable(self.variable_pres)
            model_part.AddNodalSolutionStepVariable(self.variable_trac)
            for node in self.kratos_sub_model_part_list[index].Nodes:
                model_part.CreateNewNode(node_id, node.X0, node.Y0, node.Z0)
                node_id +=1
        node_id = 1
        for index, interface_name in enumerate(output_interface_names):
            model_part = self.model.CreateModelPart(interface_name)
            model_part.AddNodalSolutionStepVariable(self.variable_disp)
            for node in self.kratos_sub_model_part_list[index].Nodes:
                model_part.CreateNewNode(node_id, node.X0, node.Y0, node.Z0)
                node_id +=1

        elem_id = 1
        for index, interface_name in enumerate(input_interface_names):
            model_part = self.model.GetModelPart(interface_name)
            for cond in self.kratos_sub_model_part_list[index].Conditions:
                node_ids = []
                for node in cond.GetNodes():
                    node_ids.append(node.Id)
                model_part.CreateNewElement('Quad',elem_id, node_ids, 0)
                elem_id +=1
        elem_id = 1
        for index, interface_name in enumerate(output_interface_names):
            model_part = self.model.GetModelPart(interface_name)
            for cond in self.kratos_sub_model_part_list[index].Conditions:
                node_ids = []
                for node in cond.GetNodes():
                    node_ids.append(node.Id)
                model_part.CreateNewElement('Quad',elem_id, node_ids, 0)
                elem_id += 1

        # # Interfaces
        self.interface_input = Interface(self.model, self.settings["interface_input"])
        self.interface_output = Interface(self.model, self.settings["interface_output"])


    def Initialize(self):
        super().Initialize()



    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.structural_analysis.time = self.structural_analysis.solver.AdvanceInTime(self.structural_analysis.time)
        self.structural_analysis.InitializeSolutionStep()
        self.structural_analysis.solver.Predict()
        self.iteration = 0




    def SolveSolutionStep(self, interface_input):

        self.interface_input.SetPythonList(interface_input.GetPythonList())
        self.ApplyLoad()
        self.structural_analysis.solver.SolveSolutionStep()
        interface_output = self.GetDisplacementInterface()
        self.iteration += 1
        return interface_output

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.structural_analysis.FinalizeSolutionStep()
        #self.structural_analysis.OutputSolutionStep()


    def Finalize(self):
        pass

    def GetInterfaceInput(self):
        return self.interface_input.deepcopy()

    def SetInterfaceInput(self):
        Exception("This solver interface provides no mapping.")

    def GetInterfaceOutput(self):
        return self.interface_output.deepcopy()

    def SetInterfaceOutput(self):
        Exception("This solver interface provides no mapping.")

    def ApplyLoad(self):
        input_interface_names = self.settings["interface_input"].keys()

        for index, interface_name in enumerate(input_interface_names):
            model_part = self.model.GetModelPart(interface_name)
            kratos_nodes = self.kratos_sub_model_part_list[index].Nodes
            node_id = 1

            for kratos_node in kratos_nodes:
                pressure = model_part.GetNode(node_id).GetSolutionStepValue(self.variable_pres)
                traction = model_part.GetNode(node_id).GetSolutionStepValue(self.variable_trac)
                kratos_node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE, 0, -1*pressure )
                kratos_node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.SURFACE_LOAD, 0, traction)
                node_id += 1



    def GetDisplacementInterface(self):
        output_interface_names = self.settings["interface_output"].keys()

        for index, interface_name in enumerate(output_interface_names):
            model_part = self.model.GetModelPart(interface_name)
            kratos_nodes = self.kratos_sub_model_part_list[index].Nodes
            node_id = 1
            for kratos_node in kratos_nodes:
                displacement_kratos = kratos_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                displacement = list(displacement_kratos)
                model_part.GetNode(node_id).SetSolutionStepValue(self.variable_disp, 0, displacement)
                node_id += 1
        return self.interface_output.deepcopy()




    def CheckInterface(self):

        input_interface_list = self.settings["interface_input"].keys()
        output_interface_list = self.settings["interface_output"].keys()

        if len(input_interface_list) == len(self.kratos_sub_model_part_list) and len(output_interface_list) == len(self.kratos_sub_model_part_list):
            pass
        else:
            raise Exception('Number of entries in the sub_model_parts_from_str_solver should match the number of input and output interfaces')







