from coconut import data_structure
from coconut import tools
from coconut.coupling_components.component import Component
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver

import time
import os
import numpy as np
from datetime import datetime
import socket

def create(parameters):
    return CoupledSolverRunSingleSolver(parameters)


class CoupledSolverRunSingleSolver(CoupledSolver):

    def __init__(self, parameters):
        """
        Initializes a single solver wrapper, but with the full
        save/restart and results-handling logic of the base CoupledSolver.
        """
        Component.__init__(self)

        self.parameters = parameters
        self.settings = parameters['settings']
        self.start_init_time = time.time()  # start of initialization

        if 'test_settings' not in self.parameters.keys():  # requires a new parameter input 'test_settings'
            raise KeyError('The coupled_solver "test_single_solver" requires "test_settings" which was not detected.')
        single_solver_settings = parameters['test_settings']

        # deprecated parameters
        if 'save_results' in self.settings:
            tools.print_info(f'WARNING: parameter "save_results" in {self.__class__.__name__} is deprecated,'
                             f' use "write_results" instead', layout='warning')
            if 'write_results' not in self.settings:
                self.settings['write_results'] = self.settings['save_results']

        # --- Read parameters specific to this class ---
        self.solver_index = single_solver_settings['solver_index']  # solver to be tested; starts at 0
        self.dummy_class = single_solver_settings.get('dummy_class')

        # --- Read parameters ---
        self.case_name = self.settings.get('case_name', 'case')  # case name
        self.settings['case_name'] = self.case_name  # make sure a case name is present
        self.timestep_start_global = self.settings.get('timestep_start', 0)
        self.timestep_start_current = self.settings.get('timestep_start', 0)
        self.restart = self.timestep_start_current != 0  # true if restart
        self.save_restart = self.settings.get('save_restart', -1)  # time step interval to save restart data
        self.settings['save_restart'] = self.save_restart  # in order to pass on default value
        self.write_results = self.settings.get('write_results', 0)  # time step interval to write coupling results
        self.anonymous = self.settings.get('anonymous', False)  # disables saving 'info' in the pickle file
        self.time_step = self.timestep_start_current  # time step
        self.delta_t = self.settings['delta_t']  # time step size
        tools.print_info(f'Using delta_t = {self.delta_t} and timestep_start = {self.timestep_start_current}')

        # create dummy components
        self.predictor = DummyComponent()
        self.convergence_criterion = DummyComponent()
        self.dummy_solver = None

        # solver wrapper settings
        parameters = self.parameters['solver_wrappers'][self.solver_index]
        if parameters['type'] == 'solver_wrappers.mapped':
            parameters = parameters['settings']['solver_wrapper']  # for mapped solver: the solver_wrapper itself tested
        settings = parameters['settings']

        # add delta_t and timestep_start to solver_wrapper settings
        tools.pass_on_parameters(self.settings, parameters['settings'],
                                 ['timestep_start', 'number_of_timesteps', 'delta_t', 'save_restart'])

        self.solver_wrapper = tools.create_instance(parameters)
        self.solver_wrappers = [self.solver_wrapper]  # used for printing summary

        self.components = [self.solver_wrapper]  # will only contain 1 solver wrapper

        # --- Initialize other variables ---
        self.x = None
        self.y = None
        self.time_step = self.timestep_start_current
        self.iteration = None  # iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)
        self.init_time = None
        self.start_run_time = None
        self.run_time = None
        self.run_time_previous = 0
        self.time_allocation = {'previous_calculations': []}
        self.iterations = []

        # --- Restart data loading ---
        if self.restart:
            self.restart_case = self.settings.get('restart_case', self.case_name)  # case to restart from
            self.restart_data = self.load_restart_data()
            self.restart_predictor = None  # indicates if predictor has to be restarted

        # save results variables
        if self.write_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.info = None
            self.case_name = self.settings.get('case_name', 'case')  # case name

        if self.settings.get('debug', False):
            tools.print_info(f'{self.__class__.__name__} has no debug mode: '
                             f'debug setting set to False', layout='warning')
        self.debug = False

    def initialize(self, print_components=True):
        Component.initialize(self)

        self.solver_wrapper.initialize()

        # initialize test_class
        interface_input = self.solver_wrapper.interface_input
        if self.dummy_class is None:
            self.dummy_solver = None
            tools.print_info('No test class specified, zero input will be used')
            for model_part_name, variable in interface_input.model_part_variable_pairs:
                if data_structure.variables_dimensions[variable] == 1:
                    tools.print_info(f'\t0 is used as {variable} input to {model_part_name}')
                elif data_structure.variables_dimensions[variable] == 3:
                    tools.print_info(f'\t[0 0 0] is used as {variable} input to {model_part_name}')
        else:
            if not os.path.isfile('dummy_solver.py'):
                raise ModuleNotFoundError(f'Test class specified, but no file named dummy_solver.py in {os.getcwd()}')
            module = tools.import_module('dummy_solver', 'dummy_solver.py')
            if not hasattr(module, self.dummy_class):
                raise NameError(f'Module dummy_solver has no class {self.dummy_class}')
            self.dummy_solver = getattr(module, self.dummy_class)()
            initialize_dummy_solver = getattr(self.dummy_solver, 'initialize', None)
            if callable(initialize_dummy_solver):
                initialize_dummy_solver(self.solver_wrapper.get_interface_input(), self.solver_index)
            tools.print_info(f'The functions from {self.dummy_class} will be used to calculate the following inputs:')
            for model_part_name, variable in interface_input.model_part_variable_pairs:
                if data_structure.variables_dimensions[variable] == 1:
                    tools.print_info(f'\t{variable} [Scalar] on {model_part_name}')
                elif data_structure.variables_dimensions[variable] == 3:
                    tools.print_info(f'\t{variable} [3D array] on {model_part_name}')
        tools.print_info()

        # initialize variables
        if self.solver_index == 1:
            self.x = self.solver_wrapper.get_interface_output()
            self.y = self.solver_wrapper.get_interface_input()
        else:
            self.x = self.solver_wrapper.get_interface_input()
            self.y = self.solver_wrapper.get_interface_output()

        # restart
        if self.restart:
            if not (self.x.has_same_model_parts(self.restart_data['interface_x']) and
                    self.y.has_same_model_parts(self.restart_data['interface_y'])):
                raise ValueError('Restart not possible because model parts changed')

        if self.write_results:
            results_data = None
            if self.restart:
                results_data = self.load_results_data()  # This loads data from _results.pickle
            if results_data is None:  # no results file to append to
                self.info = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} : ' \
                            f'start calculation of time step {self.timestep_start_current} on {socket.gethostname()}\n'
            self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
            self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)

        self.start_run_time = time.time()  # start of calculation
        self.init_time = self.start_run_time - self.start_init_time  # duration of initialization

        # Added restart header print from parent
        if self.restart:
            tools.print_info(80 * '═' + f'\n\tRestart from time step {self.timestep_start_current}\n' + 80 * '═')

    def solve_solution_step(self):
        interface_input = self.solver_wrapper.interface_input
        # generation of the input data
        if self.dummy_solver is not None:
            get_input_dummy_solver = getattr(self.dummy_solver, 'get_input', None)
            if callable(get_input_dummy_solver):
                interface_input = self.solver_wrapper.get_interface_input()
                get_input_dummy_solver(interface_input, self.time_step)
            else:
                for model_part_name, variable in interface_input.model_part_variable_pairs:
                    model_part = interface_input.get_model_part(model_part_name)
                    data = [getattr(self.dummy_solver, f'calculate_{variable}')(model_part.x0[i], model_part.y0[i],
                                                                                model_part.z0[i], self.time_step)
                            for i in range(model_part.size)]
                    interface_input.set_variable_data(model_part_name, variable, np.array(data))
        # store data in self.x and self.y
        if self.solver_index == 1:
            self.y = interface_input
            self.x = self.solver_wrapper.solve_solution_step(interface_input.copy()).copy()
        else:
            self.x = interface_input
            self.y = self.solver_wrapper.solve_solution_step(interface_input.copy()).copy()
        self.finalize_iteration(self.x * 0)

    def print_header(self):
        if self.time_step == self.timestep_start_current + 1:
            header = f'════════════════════════════════════════════════════════════════════════════════\n' \
                f"{'Time step':<16}{'Norm x':<28}{'Norm y':<28}"
            tools.print_info(header, flush=True)

    def print_iteration_info(self, r):
        info = f'{self.time_step:<16d}{self.x.norm():<28.17e}{self.y.norm():<28.17e}'
        tools.print_info(' │' * self.solver_level, info, flush=True)


class DummyComponent:
    def update(self, x):
        pass

    def initialize(self, *args, **kwargs):
        """Dummy initialize method."""
        pass

    def predict(self, x):
        """Dummy predict method."""
        return x.copy()

    def is_satisfied(self):
        """Dummy convergence check, always returns True."""
        return True

    def save_restart_data(self):
        """Dummy restart data method."""
        return None
