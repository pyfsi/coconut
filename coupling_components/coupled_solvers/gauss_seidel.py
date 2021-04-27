from coconut import tools
from coconut.tools import create_instance
from coconut.coupling_components.component import Component

import numpy as np
import time
import pickle
import os
import re


def create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']

        # read parameters
        self.timestep_start = self.settings['timestep_start']  # time step where calculation is started
        self.time_step = self.timestep_start  # time step
        self.delta_t = self.settings['delta_t']  # time step size

        self.predictor = create_instance(self.parameters['predictor'])
        self.convergence_criterion = create_instance(self.parameters['convergence_criterion'])
        self.solver_wrappers = []
        self.index_mapped = None
        self.index_other = None
        for index in range(2):
            parameters = self.parameters['solver_wrappers'][index]
            # add timestep_start and delta_t to solver_wrapper settings
            tools.pass_on_parameters(self.settings, parameters['settings'], ['timestep_start', 'delta_t'])
            self.solver_wrappers.append(create_instance(parameters))
            # determine index of mapped solver if present
            if parameters['type'] == 'solver_wrappers.mapped':
                self.index_mapped = index
            else:
                self.index_other = index
        if self.index_other is None:
            raise ValueError('Not both solvers may be mapped solvers.')

        self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = None  # input interface of solver 0
        self.y = None  # input interface of solver 1
        self.iteration = None  # iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)
        self.start_time = None
        self.elapsed_time = None
        self.iterations = []

        # restart
        self.save_restart = self.settings.get('save_restart', -1)  # time step interval to save restart data
        if self.timestep_start != 0:  # restart
            self.restart_case = self.settings['restart_case']  # case name to restart from
            self.restart_data = self.load_restart_data()

        # save results variables
        self.save_results = self.settings.get('save_results', False)  # set True in order to save for every iteration
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []

        # case name
        self.case_name = self.settings.get('name', 'case')  # case name
        listdir = os.listdir(os.getcwd())
        match = [re.search(r'(\w+)_(?:restart_ts\d+|results).pickle', f) for f in listdir]  # match file names
        present_case_names = [f.group(1) for f in match if not f is None]  # retain present case names
        if self.case_name in present_case_names:
            i = 1
            while self.case_name + str(i) in present_case_names:
                i += 1
            self.case_name += str(i)

        self.debug = False  # save results each iteration including residual interfaces
        if self.debug:
            self.complete_solution_r = None

    def initialize(self):
        super().initialize()

        # initialize mappers if required
        if self.index_mapped is not None:
            self.solver_wrappers[self.index_other].initialize()
            interface_input_from = self.solver_wrappers[self.index_other].get_interface_output()
            interface_output_to = self.solver_wrappers[self.index_other].get_interface_input()
            self.solver_wrappers[self.index_mapped].initialize(interface_input_from, interface_output_to)
        else:
            self.solver_wrappers[0].initialize()
            self.solver_wrappers[1].initialize()

        self.x = self.solver_wrappers[1].get_interface_output().copy()
        self.y = self.solver_wrappers[0].get_interface_output().copy()
        self.convergence_criterion.initialize()
        self.predictor.initialize(self.x)
        self.start_time = time.time()

        tools.print_info(80 * '═' + f'\n\t{self.case_name}\n' + 80 * '═')

        # restart
        if self.timestep_start != 0:
            if not (self.x.has_same_model_parts(self.restart_data['interface_x']) and
                    self.y.has_same_model_parts(self.restart_data['interface_y'])):
                raise ValueError('Restart not possible because model parts changed')
            self.predictor = self.restart_data['predictor']
            self.components[0] = self.predictor

        # update save results
        if self.save_results:
            if self.debug:
                self.complete_solution_x = np.empty((self.x.get_interface_data().shape[0], 0))
                self.complete_solution_y = np.empty((self.y.get_interface_data().shape[0], 0))
                self.complete_solution_r = np.empty((self.x.get_interface_data().shape[0], 0))
            elif self.timestep_start == 0:  # no restart
                self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
                self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)
            else:  # restart
                self.complete_solution_x = self.restart_data['interface_x'].get_interface_data().reshape(-1, 1)
                self.complete_solution_y = self.restart_data['interface_y'].get_interface_data().reshape(-1, 1)

    def initialize_solution_step(self):
        super().initialize_solution_step()

        # update time step and iteration
        self.time_step += 1
        self.iteration = 0

        # print time step
        if not self.solver_level:
            self.print_header()

        for component in self.components:
            component.initialize_solution_step()

        # update save results
        if self.save_results:
            self.residual.append([])

    def solve_solution_step(self):
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        y = self.solver_wrappers[0].solve_solution_step(self.x.copy())
        self.y = y.copy()
        xt = self.solver_wrappers[1].solve_solution_step(y)
        r = xt - self.x
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            self.x += r
            y = self.solver_wrappers[0].solve_solution_step(self.x)
            self.y = y.copy()
            xt = self.solver_wrappers[1].solve_solution_step(y)
            r = xt - self.x
            self.finalize_iteration(r)

    def finalize_iteration(self, r):
        self.iteration += 1  # increment iteration
        self.convergence_criterion.update(r)  # update convergence criterion
        self.print_iteration_info(r)  # print iteration information

        # update save results
        if self.save_results:
            self.residual[self.time_step - self.timestep_start - 1].append(r.norm())
            if self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x, self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y, self.y.get_interface_data().reshape(-1, 1)))
                self.complete_solution_r = np.hstack((self.complete_solution_r, r.get_interface_data().reshape(-1, 1)))
                self.output_solution_step()

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.predictor.update(self.x)
        for component in self.components:
            component.finalize_solution_step()

        # save data for restart
        if self.save_restart != 0 and self.time_step % self.save_restart == 0:
            output = {'predictor': self.predictor, 'interface_x': self.x, 'interface_y': self.y,
                      'parameters': {key: self.parameters[key] for key in ('type', 'settings', 'predictor')},
                      'delta_t': self.delta_t, 'time_step': self.time_step}
            self.add_restart_data(output)
            with open(self.case_name + f'_restart_ts{self.time_step}.pickle', 'wb') as file:
                pickle.dump(output, file)
        if self.save_restart < 0:
            try:
                os.remove(self.case_name + f'_restart_ts{self.time_step + self.save_restart}.pickle')
            except OSError:
                pass

        # update save results
        self.iterations.append(self.iteration)
        if self.save_results:
            if not self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x, self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y, self.y.get_interface_data().reshape(-1, 1)))

    def output_solution_step(self):
        super().output_solution_step()

        self.elapsed_time = time.time() - self.start_time
        if self.save_results:
            output = {'solution_x': self.complete_solution_x, 'solution_y': self.complete_solution_y,
                      'interface_x': self.x, 'interface_y': self.y, 'iterations': self.iterations,
                      'time': self.elapsed_time, 'residual': self.residual, 'delta_t': self.delta_t,
                      'timestep_start': self.timestep_start}
            if self.debug:
                output.update({'solution_r': self.complete_solution_r})
            with open(self.case_name + '_results.pickle', 'wb') as file:
                pickle.dump(output, file)

        for component in self.components:
            component.output_solution_step()

    def finalize(self):
        super().finalize()

        # print summary header
        if self.solver_level == 0:
            out = '╔' + 79 * '═' + '\n║\tSummary\n╠' + 79 * '═'
            tools.print_info(out)

        self.print_summary()

        for component in self.components:
            component.finalize()

    def load_restart_data(self):
        restart_file_name = self.restart_case + f'_restart_ts{self.timestep_start}.pickle'
        if not restart_file_name in os.listdir(os.getcwd()):
            raise FileNotFoundError(f'Not able to perform restart because {restart_file_name} '
                                    f'not found in {os.getcwd()}')
        else:
            with open(restart_file_name, 'rb') as restart_file:
                restart_data = pickle.load(restart_file)
            for key in ('predictor', 'type'):
                if self.parameters[key] != restart_data['parameters'][key]:
                    raise ValueError(f'Restart not possible because {key} changed in coupled solver')
            self.add_restart_check(restart_data)
            if self.delta_t != restart_data['delta_t']:
                tools.print_info('Time step size has changed upon restart', layout='warning')
        return restart_data

    def add_restart_check(self, restart_data):
        pass

    def add_restart_data(self, restart_data):
        pass

    def print_summary(self):
        solver_run_times = []
        pre = '║' + ' │' * self.solver_level
        out = ''
        if self.solver_level == 0:
            out += f'{pre}Elapsed time: {self.elapsed_time:0.3f}s\n'
        out += f'{pre}Percentage of total calculation time:\n'
        for solver in self.solver_wrappers:
            solver_run_times.append(solver.run_time / self.elapsed_time * 100)
            out += f'{pre}\t{solver.__class__.__name__}: {solver_run_times[-1]:0.1f}%\n'
            if solver.__class__.__name__ == 'SolverWrapperMapped':
                out += f'{pre}\t└─{solver.solver_wrapper.__class__.__name__}: ' \
                       f'{solver.solver_wrapper.run_time / self.elapsed_time * 100:0.1f}%\n'
        if self.solver_level == 0:
            out += f'{pre}\tCoupling: {100 - sum(solver_run_times):0.1f}%\n'
        out += f'{pre}Average number of iterations per time step: {np.array(self.iterations).mean():0.2f}'
        if self.solver_level == 0:
            out += '\n╚' + 79 * '═'
        tools.print_info(out)

    def print_header(self):
        header = (80 * '═' + f'\n\tTime step {self.time_step}\n' +
                  80 * '═' + f'\n{"Iteration":<16}{"Norm residual":<28}')
        tools.print_info(header, flush=True)

    def print_iteration_info(self, r):
        info = f'{self.iteration:<16d}{r.norm():<28.17e}'
        tools.print_info(' │' * self.solver_level, info, flush=True)

    def print_components_info(self, pre):
        tools.print_info(pre, 'The coupled solver ', self.__class__.__name__, ' has the following components:')
        tools.print_components_info(pre, self.components)

        # restart
        if self.timestep_start != 0:
            tools.print_info(80 * '═' + f'\n\tRestart from time step {self.timestep_start}\n' + 80 * '═')
