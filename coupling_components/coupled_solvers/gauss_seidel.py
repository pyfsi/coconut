from coconut import tools
from coconut.tools import create_instance
from coconut.coupling_components.component import Component

import numpy as np
import time
import pickle
import os
from datetime import datetime
import socket


def create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters['settings']
        self.start_init_time = time.time()  # start of initialization

        # read parameters
        self.case_name = self.settings.get('case_name', 'case')  # case name
        self.settings['case_name'] = self.case_name  # make sure a case name is present
        self.timestep_start_global = self.settings['timestep_start']  # time step for global calculation (restart)
        self.timestep_start_current = self.settings['timestep_start']  # time step start for this calculation (restart)
        self.restart = self.timestep_start_current != 0  # true if restart
        self.save_restart = self.settings.get('save_restart', -1)  # time step interval to save restart data
        self.settings['save_restart'] = self.save_restart  # in order to pass on default value
        self.save_results = self.settings.get('save_results', 0)  # time step interval to save results
        self.time_step = self.timestep_start_current  # time step
        self.delta_t = self.settings['delta_t']  # time step size

        self.predictor = create_instance(self.parameters['predictor'])
        self.convergence_criterion = create_instance(self.parameters['convergence_criterion'])
        self.solver_wrappers = []
        self.index_mapped = None
        self.index_other = None
        for index in range(2):
            parameters = self.parameters['solver_wrappers'][index]
            # add timestep_start, delta_t and save_restart to solver_wrapper settings
            tools.pass_on_parameters(self.settings, parameters['settings'], ['timestep_start', 'delta_t',
                                                                             'save_restart','number_of_timesteps'])
            self.solver_wrappers.append(create_instance(parameters))
            # determine index of mapped solver if present
            if parameters['type'] == 'solver_wrappers.mapped' or parameters['type'] == 'solver_wrappers.mapped_updated' :
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
        self.init_time = None
        self.start_run_time = None
        self.run_time = None
        self.run_time_previous = 0
        self.iterations = []

        # restart
        if self.restart:
            self.restart_case = self.settings.get('restart_case', self.case_name)  # case to restart from
            self.restart_data = self.load_restart_data()

        # save results variables
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.info = None

        self.debug = self.settings.get('debug', False)  # save results each iteration including residual interfaces
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
        self.start_run_time = time.time()  # start of calculation
        self.init_time = self.start_run_time - self.start_init_time  # duration of initialization

        if self.solver_level == 0:
            title = '╔' + 78 * '═' + f'╗\n║{self.case_name.upper():^78}║\n╚' + 78 * '═' + '╝\n'
            tools.print_info(title)

        # restart
        if self.restart:
            if not (self.x.has_same_model_parts(self.restart_data['interface_x']) and
                    self.y.has_same_model_parts(self.restart_data['interface_y'])):
                raise ValueError('Restart not possible because model parts changed')
            self.predictor = self.restart_data['predictor']
            self.components[0] = self.predictor

        # update save results
        if self.save_results:
            results_data = None
            if self.restart:
                results_data = self.load_results_data()
            if results_data is None:  # no results file to append to
                self.info = f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} : ' \
                    f'start calculation of time step {self.timestep_start_current} on {socket.gethostname()}\n'
                if self.debug:
                    self.complete_solution_x = np.empty((self.x.get_interface_data().shape[0], 0))
                    self.complete_solution_y = np.empty((self.y.get_interface_data().shape[0], 0))
                    self.complete_solution_r = np.empty((self.x.get_interface_data().shape[0], 0))
                else:
                    self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
                    self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)

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
            self.residual[self.time_step - self.timestep_start_global - 1].append(r.norm())
            if self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x,
                                                      self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y,
                                                      self.y.get_interface_data().reshape(-1, 1)))
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
            if self.save_restart < 0 and self.time_step + self.save_restart > self.timestep_start_current:
                try:
                    os.remove(self.case_name + f'_restart_ts{self.time_step + self.save_restart}.pickle')
                except OSError:
                    pass

        # update save results
        self.iterations.append(self.iteration)
        if self.save_results:
            if not self.debug:
                self.complete_solution_x = np.hstack((self.complete_solution_x,
                                                      self.x.get_interface_data().reshape(-1, 1)))
                self.complete_solution_y = np.hstack((self.complete_solution_y,
                                                      self.y.get_interface_data().reshape(-1, 1)))

    def output_solution_step(self):
        super().output_solution_step()

        self.run_time = time.time() - self.start_run_time  # duration of calculation
        if self.save_results != 0 and (self.time_step % self.save_results == 0 or
                                       (self.save_restart != 0 and self.time_step % self.save_restart == 0)):
            output = {'solution_x': self.complete_solution_x, 'solution_y': self.complete_solution_y,
                      'interface_x': self.x, 'interface_y': self.y, 'iterations': self.iterations,
                      'run_time': self.run_time + self.run_time_previous, 'residual': self.residual,
                      'delta_t': self.delta_t, 'timestep_start': self.timestep_start_global,
                      'case_name': self.case_name, 'info': self.info}
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

        for component in self.components:
            component.finalize()

        self.print_summary()

    def load_restart_data(self):
        restart_file_name = self.restart_case + f'_restart_ts{self.timestep_start_current}.pickle'
        if restart_file_name not in os.listdir(os.getcwd()):
            raise FileNotFoundError(f'Not able to perform restart because {restart_file_name} '
                                    f'not found in {os.getcwd()}')
        else:
            with open(restart_file_name, 'rb') as restart_file:
                restart_data = pickle.load(restart_file)
            for key in ('predictor', 'type'):
                if self.parameters[key] != restart_data['parameters'][key]:
                    raise ValueError(f'Restart not possible because {key} changed in coupled solver')
            self.check_restart_data(restart_data)
            if self.delta_t != restart_data['delta_t']:
                raise ValueError(f"Time step size has changed upon restart:\n\told: {restart_data['delta_t']}s"
                                 f"\n\tnew: {self.delta_t}s")
        return restart_data

    def check_restart_data(self, restart_data):
        pass

    def add_restart_data(self, restart_data):
        pass

    def load_results_data(self):
        results_file_name = f'{self.case_name}_results.pickle'
        try:
            with open(results_file_name, 'rb') as results_file:
                results_data = pickle.load(results_file)
        except FileNotFoundError:
            tools.print_info(f'Not able to append results to {results_file_name} because file not found\n'
                             f' Saving results to new file: {results_file_name}', layout='warning')
            return
        if self.debug != ('solution_r' in results_data.keys()):
            raise ValueError(f'Value of debug attribute in {self.__class__.__name__} can not be changed upon restart')
        self.timestep_start_global = results_data['timestep_start']
        self.complete_solution_x = results_data['solution_x'][:, :self.timestep_start_current
                                                              - self.timestep_start_global + 1]
        self.complete_solution_y = results_data['solution_y'][:, :self.timestep_start_current
                                                              - self.timestep_start_global + 1]
        self.x = results_data['interface_x']
        self.y = results_data['interface_y']
        self.iterations = results_data['iterations'][:self.timestep_start_current - self.timestep_start_global]
        self.run_time_previous = results_data['run_time']
        self.residual = results_data['residual'][:self.timestep_start_current - self.timestep_start_global]
        self.info = results_data.get('info', '') + '' + f'{datetime.now().strftime("%Y-%m-%d %H:%M:%S")} :' \
            f' restart calculation from time step {self.timestep_start_current} on {socket.gethostname()}\n'
        if self.debug:
            tools.print_info(f'Restart in debug mode may not append results to pickle file correctly', layout='warning')
            self.complete_solution_r = results_data['solution_r']
        return results_data

    def print_summary(self):
        solver_init_time_percs = []
        solver_run_time_percs = []
        pre = '║' + ' │' * self.solver_level
        out = ''
        if self.solver_level == 0:
            out += f'{pre}Total calculation time{" (after restart)" if self.restart else ""}:' \
                   f' {self.init_time + self.run_time:.3f}s\n'
        # initialization time
        if self.solver_level == 0:
            out += f'{pre}Initialization time: {self.init_time:0.3f}s\n'
        out += f'{pre}Distribution of initialization time:\n'
        for solver in self.solver_wrappers:
            solver_init_time_percs.append(solver.init_time / self.init_time * 100)
            out += f'{pre}\t{solver.__class__.__name__}: {solver.init_time:.0f}s ({solver_init_time_percs[-1]:0.1f}%)\n'
            if solver.__class__.__name__ == 'SolverWrapperMapped':
                out += f'{pre}\t└─{solver.solver_wrapper.__class__.__name__}: {solver.solver_wrapper.init_time:.0f}s' \
                       f' ({solver.solver_wrapper.init_time / self.init_time * 100:0.1f}%)\n'
        if self.solver_level == 0:
            out += f'{pre}\tOther: {self.init_time - sum([s.init_time for s in self.solver_wrappers]):.0f}s' \
                   f' ({100 - sum(solver_init_time_percs):0.1f}%)\n'
        # run time
        if self.solver_level == 0:
            out += f'{pre}Run time{" (after restart)" if self.restart else ""}: {self.run_time:0.3f}s\n'
        out += f'{pre}Distribution of run time:\n'
        for solver in self.solver_wrappers:
            solver_run_time_percs.append(solver.run_time / self.run_time * 100)
            out += f'{pre}\t{solver.__class__.__name__}: {solver.run_time:.0f}s ({solver_run_time_percs[-1]:0.1f}%)\n'
            if solver.__class__.__name__ == 'SolverWrapperMapped':
                out += f'{pre}\t└─{solver.solver_wrapper.__class__.__name__}: {solver.solver_wrapper.run_time:.0f}s' \
                       f' ({solver.solver_wrapper.run_time / self.run_time * 100:0.1f}%)\n'
        if self.solver_level == 0:
            out += f'{pre}\tCoupling: {self.run_time - sum([s.run_time for s in self.solver_wrappers]):.0f}s' \
                   f' ({100 - sum(solver_run_time_percs):0.1f}%)\n'
        out += f'{pre}Average number of iterations per time step' \
               f'{" (including before restart)" if self.restart else ""}: {np.array(self.iterations).mean():0.2f}'
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
        if self.restart:
            tools.print_info(80 * '═' + f'\n\tRestart from time step {self.timestep_start_current}\n' + 80 * '═')
