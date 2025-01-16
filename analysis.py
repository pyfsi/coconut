from coconut.tools import create_instance, pass_on_parameters
import json
import sys


class Analysis:
    def __init__(self, parameters):
        self.parameters = parameters
        self.settings = parameters['settings']

        self.number_of_timesteps = self.settings['number_of_timesteps']
        self.settings['save_restart'] = self.settings.get('save_restart', -1)  # restart data save interval, default -1
        pass_on_parameters(self.settings, parameters['coupled_solver']['settings'],
                           ['timestep_start', 'number_of_timesteps', 'delta_t', 'save_restart'])

        # Python version: 3.6 or higher
        if sys.version_info < (3, 6):
            raise RuntimeError('Python version 3.6 or higher required.')

        self.coupled_solver = create_instance(self.parameters['coupled_solver'])

    def run(self):
        self.initialize()
        self.run_solution_loop()
        self.finalize()

    def initialize(self):
        self.coupled_solver.initialize()

    def run_solution_loop(self):
        for _ in range(self.number_of_timesteps):
            self.coupled_solver.initialize_solution_step()
            self.coupled_solver.solve_solution_step()
            self.coupled_solver.finalize_solution_step()
            self.coupled_solver.output_solution_step()

    def finalize(self):
        self.coupled_solver.finalize()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='The main file which runs the coupled simulation')
    parser.add_argument('parameter_file_name', help='Name of the parameter file containing coupling information',
                        type=str, nargs=1)

    # import parameters using the data structure
    with open(parser.parse_args().parameter_file_name, 'r') as parameter_file:
        parameters = json.load(parameter_file)

    simulation = Analysis(parameters)
    simulation.run()
