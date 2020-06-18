from coconut import data_structure
from coconut.coupling_components.tools import CreateInstance


class Analysis:
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]

        self.number_of_timesteps = self.settings["number_of_timesteps"].GetInt()

        self.coupled_solver = CreateInstance(self.parameters["coupled_solver"])

    def Run(self):
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def Initialize(self):
        self.coupled_solver.Initialize()
        self.coupled_solver.Check()
        self.coupled_solver.PrintInfo(' ')

    def RunSolutionLoop(self):
        for _ in range(self.number_of_timesteps):
            self.coupled_solver.InitializeSolutionStep()
            self.coupled_solver.SolveSolutionStep()
            self.coupled_solver.FinalizeSolutionStep()
            self.coupled_solver.OutputSolutionStep()

    def Finalize(self):
        self.coupled_solver.Finalize()


if __name__ == '__main__':
    from sys import argv

    # Check number of command line arguments
    if len(argv) != 2:
        err_msg = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python analysis.py <parameter>.json"\n'
        raise Exception(err_msg)

    # Import data structure
    parameter_file_name = argv[1]

    # Import parameters using the data structure
    with open(parameter_file_name, 'r') as parameter_file:
        parameters = data_structure.Parameters(parameter_file.read())

    simulation = Analysis(parameters)
    simulation.Run()
