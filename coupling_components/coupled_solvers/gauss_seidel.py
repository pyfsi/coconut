from coconut import tools
from coconut.tools import create_instance
from coconut.coupling_components.component import Component

import numpy as np
import time
import pickle


def create(parameters):
    return CoupledSolverGaussSeidel(parameters)


class CoupledSolverGaussSeidel(Component):
    def __init__(self, parameters):
        super().__init__()

        self.parameters = parameters
        self.settings = parameters["settings"]

        self.timestep_start = self.settings["timestep_start"]  # time step where calculation is started
        self.n = self.timestep_start  # time step
        self.delta_t = self.settings["delta_t"]  # time step size

        self.predictor = create_instance(self.parameters["predictor"])
        self.convergence_criterion = create_instance(self.parameters["convergence_criterion"])
        self.solver_wrappers = []
        for index in range(2):
            # add timestep_start and delta_t to solver_wrapper settings
            parameters = self.parameters["solver_wrappers"][index]
            if parameters["type"] == "solver_wrappers.mapped":
                settings = parameters["settings"]["solver_wrapper"]["settings"]
            else:
                settings = parameters["settings"]

            for key in ["timestep_start", "delta_t"]:
                if key in settings:
                    tools.print_info(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
                settings[key] = self.settings[key]

            self.solver_wrappers.append(create_instance(parameters))

        self.components = [self.predictor, self.convergence_criterion, self.solver_wrappers[0], self.solver_wrappers[1]]

        self.x = None
        self.y = None
        self.iteration = None  # iteration
        self.solver_level = 0  # 0 is main solver (time step is printed)

        self.start_time = None
        self.elapsed_time = None
        self.iterations = []
        self.save_results = self.settings.get("save_results", False)  # set True in order to save for every iteration
        if self.save_results:
            self.complete_solution_x = None
            self.complete_solution_y = None
            self.residual = []
            self.case_name = self.settings.get("name", "results")  # case name

    def initialize(self):
        super().initialize()

        for component in self.components[1:]:
            component.initialize()

        # construct mappers if required
        index_mapped = None
        index_other = None
        for i in range(2):
            type = self.parameters["solver_wrappers"][i]["type"]
            if type == "solver_wrappers.mapped":
                index_mapped = i
            else:
                index_other = i
        if index_other is None:
            raise ValueError("Not both solvers may be mapped solvers.")
        if index_mapped is not None:
            # construct input mapper
            interface_input_from = self.solver_wrappers[index_other].get_interface_output()
            self.solver_wrappers[index_mapped].set_interface_input(interface_input_from)

            # construct output mapper
            interface_output_to = self.solver_wrappers[index_other].get_interface_input()
            self.solver_wrappers[index_mapped].set_interface_output(interface_output_to)

        # initialize variables
        self.x = self.solver_wrappers[1].get_interface_output()
        self.y = self.solver_wrappers[0].get_interface_output()
        self.predictor.initialize(self.x)

        if self.save_results:
            self.complete_solution_x = self.x.get_interface_data().reshape(-1, 1)
            self.complete_solution_y = self.y.get_interface_data().reshape(-1, 1)
        self.start_time = time.time()

    def initialize_solution_step(self):
        super().initialize_solution_step()

        for component in self.components:
            component.initialize_solution_step()

        self.n += 1  # increment time step
        self.iteration = 0

        # print time step
        if not self.solver_level:
            self.print_header()

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
        self.iteration += 1
        self.convergence_criterion.update(r)
        # print iteration information
        self.print_iteration_info(r)

        if self.save_results:
            self.residual[self.n - 1].append(np.linalg.norm(r.get_interface_data()))

    def finalize_solution_step(self):
        super().finalize_solution_step()

        self.iterations.append(self.iteration)
        if self.save_results:
            self.complete_solution_x = np.hstack((self.complete_solution_x, self.x.get_interface_data().reshape(-1, 1)))
            self.complete_solution_y = np.hstack((self.complete_solution_y, self.y.get_interface_data().reshape(-1, 1)))

        self.predictor.update(self.x)
        for component in self.components:
            component.finalize_solution_step()

    def output_solution_step(self):
        super().output_solution_step()

        for component in self.components:
            component.output_solution_step()

    def finalize(self):
        super().finalize()

        if self.solver_level == 0:
            out = f"╔═══════════════════════════════════════════════════════════════════════════════\n" \
                  f"║\tSummary\n" \
                  f"╠═══════════════════════════════════════════════════════════════════════════════"
            tools.print_info(out)

        for component in self.components:
            component.finalize()

        self.elapsed_time = time.time() - self.start_time
        self.print_summary()
        if self.save_results:
            output = {"solution_x": self.complete_solution_x, "solution_y": self.complete_solution_y,
                      "interface_x": self.x, "interface_y": self.y, "iterations": self.iterations,
                      "time": self.elapsed_time, "residual": self.residual, "delta_t": self.delta_t,
                      "timestep_start": self.timestep_start}
            pickle.dump(output, open(self.case_name + '.pickle', 'wb'))

    def print_summary(self):
        solver_run_times = []
        pre = '║' + ' │' * self.solver_level
        out = ""
        if self.solver_level == 0:
            out += f"{pre}Elapsed time: {self.elapsed_time:0.3f}s\n"
        out += f"{pre}Percentage of total calculation time:\n"
        for solver in self.solver_wrappers:
            solver_run_times.append(solver.run_time / self.elapsed_time * 100)
            out += f"{pre}\t{solver.__class__.__name__}: {solver_run_times[-1]:0.1f}%\n"
            if solver.__class__.__name__ == "SolverWrapperMapped":
                out += f"{pre}\t└─{solver.solver_wrapper.__class__.__name__}: " \
                       f"{solver.solver_wrapper.run_time / self.elapsed_time * 100:0.1f}%\n"
        if self.solver_level == 0:
            out += f"{pre}\tCoupling: {100 - sum(solver_run_times):0.1f}%\n"
        out += f"{pre}Average number of iterations per time step: {np.array(self.iterations).mean():0.2f}"
        if self.solver_level == 0:
            out += f"\n╚═══════════════════════════════════════════════════════════════════════════════"
        tools.print_info(out)

    def print_header(self):
        header = f"════════════════════════════════════════════════════════════════════════════════\n" \
                 f"\tTime step {self.n}\n" \
                 f"════════════════════════════════════════════════════════════════════════════════\n" \
                 f"{'Iteration':<16}{'Norm residual':<28}"
        tools.print_info(header, flush=True)

    def print_iteration_info(self, r):
        info = f"{self.iteration:<16d}{r.norm():<28.17e}"
        tools.print_info(' │' * self.solver_level, info, flush=True)

    def check(self):
        super().check()

        for component in self.components:
            component.check()

    def print_components_info(self, pre):
        tools.print_info(pre, "The coupled solver ", self.__class__.__name__, " has the following components:")
        tools.print_components_info(pre, self.components)
