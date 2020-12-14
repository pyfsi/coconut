from coconut.coupling_components.tools import create_instance
from coconut.coupling_components.component import Component
from coconut.coupling_components import tools

""" proposed changes to mapped.py
- do initialization of mappers in Initialize method, would be more logical
- remove all SetInterfaceInput/Output methods?
- use copy in GetInterfaceInput/Output methods?
    and just refer to actual solver wrapper in SolverWrapperMapped
- all Interfaces are stored in this mapper, e.g. self.interface_output_to and 3 others;
    I see no reason for this; furthermore, it is only useful to store it if you take copies all the time
- OutputSolutionStep is barely used; what's the deal with it??
"""


def Create(parameters):
    return SolverWrapperMapped(parameters)


class SolverWrapperMapped(Component):
    def __init__(self, parameters):
        super().__init__()

        # Read parameters
        self.parameters = parameters
        self.settings = parameters["settings"]

        # Create solver
        self.solver_wrapper = create_instance(self.settings["solver_wrapper"])

        # run time
        self.run_time = 0.0

    def Initialize(self):
        super().Initialize()

        self.solver_wrapper.Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.solver_wrapper.InitializeSolutionStep()

    @tools.TimeSolveSolutionStep
    def SolveSolutionStep(self, interface_input_from):
        self.interface_input_from = interface_input_from
        self.mapper_interface_input(self.interface_input_from, self.interface_input_to)
        self.interface_output_from = self.solver_wrapper.SolveSolutionStep(self.interface_input_to)
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.solver_wrapper.FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()

        self.solver_wrapper.Finalize()
        self.mapper_interface_input.Finalize()
        self.mapper_interface_output.Finalize()

    def OutputSolutionStep(self):
        super().OutputSolutionStep()

        self.solver_wrapper.OutputSolutionStep()
        self.mapper_interface_input.OutputSolutionStep()
        self.mapper_interface_output.OutputSolutionStep()

    def GetInterfaceInput(self):
        # Does not contain most recent data
        # *** shouldn't this just call the underlying solver wrapper?
        return self.interface_input_from

    def SetInterfaceInput(self, interface_input_from):
        # Create input mapper
        self.interface_input_from = interface_input_from.deepcopy()
        self.interface_input_to = self.solver_wrapper.GetInterfaceInput()

        self.mapper_interface_input = create_instance(self.settings["mapper_interface_input"])
        self.mapper_interface_input.Initialize(self.interface_input_from, self.interface_input_to)

    def GetInterfaceOutput(self):
        self.interface_output_from = self.solver_wrapper.GetInterfaceOutput()
        self.mapper_interface_output(self.interface_output_from, self.interface_output_to)
        return self.interface_output_to.deepcopy()

    def SetInterfaceOutput(self, interface_output_to):
        # Create output mapper
        self.interface_output_to = interface_output_to.deepcopy()
        self.interface_output_from = self.solver_wrapper.GetInterfaceOutput()

        self.mapper_interface_output = create_instance(self.settings["mapper_interface_output"])
        self.mapper_interface_output.Initialize(self.interface_output_from, self.interface_output_to)

    def PrintComponentsInfo(self, pre):
        tools.Print(pre, "The component ", self.__class__.__name__, " maps the following solver wrapper:")
        pre = tools.UpdatePre(pre)
        self.solver_wrapper.PrintComponentsInfo(pre + '├─')
        tools.Print(pre, '├─', "Input mapper:")
        self.mapper_interface_input.PrintComponentsInfo(pre + '│ └─')
        tools.Print(pre, '└─', "Output mapper:")
        self.mapper_interface_output.PrintComponentsInfo(pre + '  └─')
