from coconut import data_structure
from coconut.coupling_components.tools import create_instance
import json


class Analysis:
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters["settings"]

        self.number_of_timesteps = self.settings["number_of_timesteps"]

        self.coupled_solver = create_instance(self.parameters["coupled_solver"])

    def Run(self):
        self.Initialize()
        self.RunSolutionStep()
        self.Finalize()

    def Initialize(self):
        self.coupled_solver.Initialize()
        self.coupled_solver.Check()
        self.coupled_solver.PrintInfo(' ')

    def RunSolutionLoop(self):
        for _ in range(self.number_of_timesteps):
            self.coupled_solver.initializeSolutionStep()
            self.coupled_solver.SolveSolutionStep()
            self.coupled_solver.FinalizeSolutionStep()
            self.coupled_solver.OutputSolutionStep()

    def Finalize(self):
        self.coupled_solver.Finalize()



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='The main file which runs the coupled simulation')
    parser.add_argument('parameter_file_name', help='Name of the parameter file containing coupling information.', type=str, nargs=1)

    # # Check number of command line arguments
    # if len(argv) != 2:
    #     err_msg = 'Wrong number of input arguments!\n'
    #     err_msg += 'Use this script in the following way:\n'
    #     err_msg += '    "python analysis.py <parameter>.json"\n'
    #     raise Exception(err_msg)

    # Import data structure
    #parameter_file_name = argv[1]

    # Import parameters using the data structure
    with open(parser.parse_args().parameter_file_name, 'r') as parameter_file:
        parameters = json.load(parameter_file)

    simulation = Analysis(parameters)
    simulation.Run()
