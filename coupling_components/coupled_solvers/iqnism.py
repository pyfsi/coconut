from coconut import tools
from coconut.coupling_components.coupled_solvers.coupled_solver import CoupledSolver

import numpy as np
import time


def create(parameters):
    return CoupledSolverIQNISM(parameters)


class CoupledSolverIQNISM(CoupledSolver):
    def __init__(self, parameters):
        super().__init__(parameters)

        # add timestep_start, delta_t, save_restart and case_name to surrogate settings
        if 'settings' not in self.settings['surrogate']:
            self.settings['surrogate']['settings'] = {}
        tools.pass_on_parameters(self.settings, self.settings['surrogate']['settings'], ['timestep_start', 'delta_t',
                                                                                         'save_restart', 'case_name'])
        # add timestep_start, delta_t, save_restart to model settings
        tools.pass_on_parameters(self.settings, self.settings['model']['settings'],
                                 ('timestep_start', 'delta_t', 'save_restart'))

        self.model = tools.create_instance(self.settings['model'])
        self.surrogate_mapped = (self.settings['surrogate']['type'] == 'coupled_solvers.models.mapped')
        self.surrogate_dummy = (self.settings['surrogate']['type'] == 'coupled_solvers.models.dummy_model')
        self.surrogate = tools.create_instance(self.settings['surrogate'], 'models.dummy_model')
        self.omega = self.settings.get('omega', 1)  # relaxation factor
        self.surrogate_modes = self.settings.get('surrogate_modes')  # number of surrogates modes to use
        self.surrogate_synchronize = self.settings.get('surrogate_synchronize', True)  # synchronize surrogate

        if ((self.surrogate_modes == 0 or self.settings['surrogate']['type'] == 'coupled_solvers.models.dummy_model')
                and 'omega' not in self.settings):
            raise ValueError('A relaxation factor (omega) is required when no surrogate Jacobian is used')

        self.restart_model_surrogate = [False, False]  # indicates if model has to be restarted

    def initialize(self):
        super().initialize()

        self.surrogate.solver_level = self.solver_level + 1

        # initialize mapper for surrogate if required
        if self.surrogate_mapped:
            interface_input_from = self.solver_wrappers[0].get_interface_input()
            interface_output_to = self.solver_wrappers[1].get_interface_output()  # same interface a input_from
            self.surrogate.initialize(interface_input_from, interface_output_to)
        else:
            self.surrogate.initialize()

        self.model.size_in = self.model.size_out = self.x.size
        self.model.out = self.x.copy()
        self.model.initialize()

        restart_model, restart_surrogate = self.restart_model_surrogate
        if restart_model:  # restart with the same model type
            self.model.restart(self.restart_data['model'])
        if restart_surrogate:  # restart with the same surrogate type
            self.surrogate.restart(self.restart_data['surrogate'])

        self.components += [self.model, self.surrogate]

        # set initial surrogate value in surrogate predictor
        if self.surrogate.provides_get_solution and not self.restart:
            if hasattr(self.predictor, 'update_surrogate') and not self.surrogate_dummy:
                self.predictor.update_surrogate(self.surrogate.get_interface_output())

        self.start_run_time = time.time()  # reset start of calculation
        self.init_time = self.start_run_time - self.start_init_time  # reset duration of initialization

    def solve_solution_step(self):
        # solve surrogate
        if self.surrogate.provides_get_solution:
            x_surrogate = self.surrogate.get_solution()
            if hasattr(self.predictor, 'update_surrogate') and not self.surrogate_dummy:
                self.predictor.update_surrogate(x_surrogate)
        # initial value
        self.x = self.predictor.predict(self.x)
        # first coupling iteration
        self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
        xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
        r = xt - self.x
        self.model.add(r.copy(), xt.copy())
        self.surrogate.add(r.copy(), xt)  # only used when derivative info of surrogate is updated every iteration
        self.finalize_iteration(r)
        # coupling iteration loop
        while not self.convergence_criterion.is_satisfied():
            dr = -1 * r
            if not self.model.is_ready():
                if not self.surrogate.is_ready():
                    dx = -self.omega * dr
                else:
                    dx = self.surrogate.predict(dr.copy(), modes=self.surrogate_modes) - dr
                    # relax other modes
                    dx -= (self.omega - 1.0) * self.surrogate.filter_q(dr.copy(), modes=self.surrogate_modes)
            else:
                if not self.surrogate.is_ready():
                    dx = self.model.predict(dr.copy()) - dr
                else:
                    dx = self.model.predict(dr.copy()) + self.surrogate.predict(self.model.filter_q(dr.copy()),
                                                                                modes=self.surrogate_modes) - dr
            self.x += dx
            self.y = self.solver_wrappers[0].solve_solution_step(self.x.copy()).copy()
            xt = self.solver_wrappers[1].solve_solution_step(self.y.copy()).copy()
            r = xt - self.x
            self.model.add(r.copy(), xt.copy())
            self.surrogate.add(r.copy(), xt)  # only used when derivative information of surrogate is function of x
            self.finalize_iteration(r)
        # synchronize
        if self.surrogate_synchronize and self.surrogate.provides_set_solution:
            self.surrogate.set_solution(self.x.copy())

    def check_restart_data(self, restart_data, coupled_solver_settings=None):
        continue_check = super().check_restart_data(restart_data, ['omega', 'surrogate_modes', 'surrogate_synchronize'])
        if continue_check:
            for i, model_name in enumerate(['model', 'surrogate']):
                old_model_type = restart_data['parameters']['settings'][model_name]['type']
                new_model_type = self.parameters['settings'][model_name]['type']
                if new_model_type != old_model_type:
                    tools.print_info(f'Model type changed from "{old_model_type}" to "{new_model_type}"', layout='blue')
                else:
                    self.restart_model_surrogate[i] = True
            if self.restart_model_surrogate[0]:
                self.model.check_restart_data(restart_data['parameters']['settings']['model'])

    def add_restart_data(self, restart_data):
        restart_data.update({'model': self.model.save_restart_data(), 'surrogate': self.surrogate.save_restart_data()})

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
        surrogate_init_time_perc = self.surrogate.init_time / self.init_time * 100
        out += f'{pre}\t{self.surrogate.__class__.__name__}: {self.surrogate.init_time:.0f}s' \
               f' ({surrogate_init_time_perc:0.1f}%)\t\n'
        if self.surrogate.__class__.__name__ == 'ModelMapped':
            out += f'{pre}\t└─{self.surrogate.surrogate.__class__.__name__}: ' \
                   f'{self.surrogate.surrogate.init_time:.0f}s' \
                   f' ({self.surrogate.surrogate.init_time / self.init_time * 100:0.1f}%)\n'
        if self.solver_level == 0:
            other_time = self.init_time - sum([s.init_time for s in self.solver_wrappers]) \
                         - self.surrogate.init_time
            out += f'{pre}\tOther: {other_time:.0f}s' \
                   f' ({100 - sum(solver_init_time_percs) - surrogate_init_time_perc:0.1f}%)\n'
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
        surrogate_run_time_perc = self.surrogate.run_time / self.run_time * 100
        out += f'{pre}\t{self.surrogate.__class__.__name__}: {self.surrogate.run_time:.0f}s' \
               f' ({surrogate_run_time_perc:0.1f}%)\t\n'
        if self.surrogate.__class__.__name__ == 'ModelMapped':
            out += f'{pre}\t└─{self.surrogate.surrogate.__class__.__name__}: {self.surrogate.surrogate.run_time:.0f}s' \
                   f' ({self.surrogate.surrogate.run_time / self.run_time * 100:0.1f}%)\n'
        if self.solver_level == 0:
            coupling_time = self.run_time - sum([s.run_time for s in self.solver_wrappers]) \
                            - self.surrogate.run_time
            out += f'{pre}\tCoupling: {coupling_time:.0f}s' \
                   f' ({100 - sum(solver_run_time_percs) - surrogate_run_time_perc:0.1f}%)\n'
        out += f'{pre}Average number of iterations per time step' \
               f'{" (including before restart)" if self.restart else ""}: {np.array(self.iterations).mean():0.2f}'
        if self.solver_level == 0:
            out += '\n╚' + 79 * '═'
        tools.print_info(out)

    def get_time_allocation(self):
        time_allocation = super().get_time_allocation()
        for time_type in ('init_time', 'run_time'):
            time_allocation_sub = time_allocation[time_type]
            time_allocation_sub.update({f'surrogate': self.surrogate.get_time_allocation()[time_type]})
            time_allocation_sub['other'] = self.__getattribute__(time_type) - sum([
                s.__getattribute__(time_type) for s in self.solver_wrappers]) - self.surrogate.__getattribute__(
                time_type)
        return time_allocation
