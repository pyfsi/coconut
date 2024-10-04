import collections.abc
import pickle
import warnings
from fractions import Fraction

import matplotlib.animation as ani
import matplotlib.pyplot as plt
import numpy as np
from coconut import tools
from coconut.data_structure import variables_dimensions


class PostProcess:
    def __init__(self, pickle_file_path: str):
        self.pickle_file_path = pickle_file_path
        with open(self.pickle_file_path, 'rb') as f:
            self.data = pickle.load(f)
        self.debug_run = False
        self.solutions = [self.data['solution_x'], self.data['solution_y']]
        if 'solution_r' in self.data:
            self.solutions.append(self.data['solution_r'])
            self.debug_run = True

        self.dt = self.data['delta_t']
        self.time_step_start = self.data['timestep_start']
        self.time_steps = self.solutions[0].shape[1] - 1

        self.case_name = self.data['case_name']

        self.interface_names = ['interface_x', 'interface_y']
        self.interfaces = [self.data[interface] for interface in self.interface_names]

        self.subsets = []

    def print_model_parts_variables(self):
        tools.print_info(self.get_model_parts_variables())

    def get_model_parts_variables(self, pre='', print_size=False):
        # returns a string with an overview of model parts and variables
        out = ''
        for interface_name in self.interface_names:
            out += f'{pre}interface "{interface_name}" consisting of\n'
            interface = self.data[interface_name]
            for mp_dict in interface.parameters:
                mp_name = mp_dict['model_part']
                mp_size = f'of size {interface.get_model_part(mp_name).size} ' if print_size else ''
                out += f'{pre}{" " * 4}model part "{mp_name}" {mp_size}with variables\n'
                for var in mp_dict['variables']:
                    dim = variables_dimensions[var]
                    out += f'{pre}{" " * 8}{var} ({"scalar" if dim == 1 else "vector"})\n'
        return out

    def get_model_part_names(self, interface_name):
        if interface_name not in self.interface_names:
            raise ValueError(f'Interface "{interface_name}" not found\nOptions:\n{" ,".join(self.interface_names)}')
        return [mp_dict['model_part'] for mp_dict in
                self.interfaces[self.interface_names.index(interface_name)].parameters]

    def _get_tree(self, dc, total, pre=''):
        # returns a string providing an overview of a time allocation dictionary
        out = ''
        for key, value in dc.items():
            last = key == tuple(dc.keys())[-1]
            connector = '└─' if last else '├─'
            if key == 'total':
                continue
            if isinstance(value, float):
                out += f'{pre}{connector}{key.replace("_", " ").capitalize()}: {value:.0f}s ' \
                       f'({value / total * 100:0.1f}%)\n'
            elif isinstance(value, dict):
                out += f'{pre}{connector}{key.replace("_", " ").capitalize()}: {value["total"]:.0f}s ' \
                       f'({value["total"] / total * 100:0.1f}%)\n'
                out += self._get_tree(value, total, pre=pre + ('  ' if last else '│ '))
            else:
                raise ValueError(f'Unknown value type "{type(value)}"')
        return out

    def get_summary(self, pre=''):
        time_allocation = self.data['time_allocation']
        previous_calculations = time_allocation['previous_calculations']
        restart = len(previous_calculations) != 0
        out = f'{pre}Summary\n'
        out += f'{pre}Total calculation time{" (after restart)" if restart else ""}: ' \
               f'{time_allocation["total"]:.3f}s\n'

        # initialization time
        total = time_allocation["init_time"]["total"]
        out += f'{pre}Initialization time: {total:0.3f}s\n'
        out += f'{pre}Distribution of initialization time:\n'
        out += self._get_tree(time_allocation['init_time'], total, pre=pre + '  ')

        # run time and save time (part of run time)
        total = time_allocation["run_time"]["total"]
        out += f'{pre}Run time{" (after restart)" if restart else ""}: {total:0.3f}s\n'
        out += f'{pre}Distribution of run time:\n'
        out += self._get_tree(time_allocation['run_time'], total, pre=pre + '  ')

        out += f'{pre}Average number of iterations per time step' \
               f'{" (including before restart)" if restart else ""}: {np.array(self.data["iterations"]).mean():0.2f}'
        return out

    def print_summary(self):
        tools.print_info(self.get_summary())

    def print_info(self):
        tools.print_info(self.data['info'])

    def get_data(self):
        return dict(self.data)  # return copy

    def get_residual_history(self):
        return self.data['residual'].copy()  # return copy

    def get_coupling_iterations(self):
        return self.data['iterations'].copy()  # return copy

    def add_subset(self, **kwargs):
        self.subsets.append(SubSet(self, **kwargs))
        return self.subsets[-1]

    def get_subsets(self):
        return self.subsets

    def __repr__(self):
        return f'PostProcess of pickle file {self.pickle_file_path}\nContaining {self.time_steps} ' \
               f'time step{"" if self.time_steps == 1 else "s"} with size {self.dt}s\n' \
               f'Coupling\n{self.get_model_parts_variables(" " * 4, print_size=True)}'


class SubSet:
    def __init__(self, run, interface=None, model_part=None, variable=None, component=None, sort=(0, 1, 2)):
        self.run = run

        # create mapping from model part names to interface names
        mp_names = {}
        for interface_name in self.run.interface_names:
            for model_part_dict in self.run.data[interface_name].parameters:
                mp_name = model_part_dict['model_part']
                if mp_name not in mp_names:
                    mp_names.update({mp_name: [interface_name]})
                else:
                    mp_names[mp_name].append(interface_name)

        # check if correct interface and model parameter names were provided
        options = self.run.get_model_parts_variables(" " * 4)
        if interface is None:
            if model_part is None:
                if len(self.run.interfaces) != 1 or len(mp_names.keys()) != 1:
                    raise ValueError(f'Specify interface and/or model part\nOptions:\n{options}')
                else:
                    self.interface_name = self.run.interface_names[0]
                    self.model_part_name = list(mp_names.keys())[0]
            else:
                if model_part not in mp_names.keys():
                    raise ValueError(f'Model part "{model_part}" not found\nOptions:\n{options}')
                elif len(mp_names[model_part]) != 1:
                    raise ValueError(f'Model part "{model_part}" present in multiple interfaces: specify interface\n'
                                     f'Options:\n{options}')
                else:
                    self.model_part_name = model_part
                    self.interface_name = mp_names[model_part][0]
            self.interface = self.run.data[self.interface_name]
        else:
            if interface not in self.run.interface_names:
                raise ValueError(f'Interface "{interface}" not found\nOptions:\n{options}')
            else:
                self.interface_name = interface
                self.interface = self.run.data[self.interface_name]
                if model_part is None:
                    if len(set(_[0] for _ in self.interface.model_part_variable_pairs)) > 1:
                        raise ValueError(f'Specify model part in interface "{interface}"\nOptions:\n{options}')
                    else:
                        self.model_part_name = self.interface.model_part_variable_pairs[0][0]
                else:
                    if model_part not in [_[0] for _ in self.interface.model_part_variable_pairs]:
                        raise ValueError(f'Model part "{model_part}" not found in interface "{interface}"\n'
                                         f'Options:\n{options}')
                    else:
                        self.model_part_name = model_part
        self.model_part = self.interface.get_model_part(self.model_part_name)

        self.mp_var_pairs = self.interface.model_part_variable_pairs

        # get solution data NOTE TRANSPOSE: 0-axis (rows) is time axis
        if self.run.debug_run:
            self.solution_data = self.run.data['solution_r'].T
        else:
            self.solution_data = self.run.solutions[self.run.interface_names.index(self.interface_name)].T

        # get time steps
        self.dt = self.run.dt
        self.complete_num_steps = self.solution_data.shape[0] - 1  # number of times steps
        self.steps = np.arange(self.run.time_step_start, (self.run.time_step_start + self.complete_num_steps + 1))
        self.times = self.steps * self.dt

        # get coordinates
        self.initial_coordinates = np.zeros((self.model_part.size, 3))
        for j, direction in enumerate(['x0', 'y0', 'z0']):
            self.initial_coordinates[:, j] = getattr(self.model_part, direction)
        self.complete_size = self.model_part.size  # number of points

        # sort points
        if sort:
            if len(sort) != 3:
                raise ValueError(
                    'Provide the sort order as a list or tuple containing 0, 1 and 2 referring to x, y and z')
            self.argsort = np.lexsort((self.initial_coordinates[:, sort[2]],
                                       self.initial_coordinates[:, sort[1]],
                                       self.initial_coordinates[:, sort[0]]
                                       ))
            self.initial_coordinates = self.initial_coordinates[self.argsort]
        else:
            self.argsort = None

        # set variable if possible
        self.variable = None
        self.dimension = None
        self.component = None
        self.available_vars = None
        self.complete_solution = None
        self._set_variable(variable, component=component)

        # reset masks and get solution
        self.p_mask = None
        self.size = None
        self.t_mask = None
        self.num_steps = None
        self.solution = None
        self.reset_points_selection()
        self.reset_times_selection()
        self._get_solution()

    def _set_variable(self, variable, component=None):
        # check if correct variable is provided
        self.available_vars = [pair[1] for pair in self.mp_var_pairs if pair[0] == self.model_part_name]
        if 'displacement' in self.available_vars:
            self.available_vars.append('coordinates')
        if variable is None:
            if len(self.available_vars) == 1:
                self.variable = self.available_vars[0]
            elif self.variable is None:  # no variable has been specified
                return
        else:
            if variable not in self.available_vars:
                raise ValueError(f'Variable "{variable}" is not found, choose from {", ".join(self.available_vars)}')
            else:
                self.variable = variable

        # check if component is correct
        if component is not None:
            if component in ('x', 'y', 'z'):
                self.component = ('x', 'y', 'z').index(component)
            elif component in (0, 1, 2):
                self.component = component
            else:
                raise ValueError(f'Component "{component}" invalid, choose from 0, 1, 2 or "x", "y", "z"')
        else:
            self.component = None

        # find location of data
        index = 0
        start_index = end_index = None
        for mp_name, var in self.mp_var_pairs:
            mp = self.interface.get_model_part(mp_name)
            dimension = variables_dimensions[var]
            var_search = self.variable if self.variable != 'coordinates' else 'displacement'
            if var == var_search and mp_name == self.model_part_name:  # correct location
                start_index = index
                self.dimension = dimension
                end_index = index + self.dimension * self.complete_size
                if self.complete_size != mp.size:
                    raise RuntimeError('Number of points does not match')
            index += dimension * mp.size
        if index != self.solution_data.shape[1]:
            raise RuntimeError('Size of provided solution data does not match interface')
        if start_index is None:
            raise RuntimeError('Combination of variable and model part not found')

        # data structured as time steps x points x components
        self.complete_solution = self.solution_data[:, start_index: end_index].reshape(
            self.complete_num_steps + 1, self.complete_size, self.dimension)
        if self.argsort is not None:
            self.complete_solution = self.complete_solution[:, self.argsort]

    def get_available_variables(self):
        return self.available_vars

    def get_all_initial_coordinates(self):
        return np.copy(self.initial_coordinates)

    def get_initial_coordinates_selection(self):
        return np.copy(self.initial_coordinates[self.p_mask])

    def get_size(self):
        return self.size

    def get_all_times(self):
        return np.copy(self.times)

    def get_times_selection(self):
        return np.copy(self.times[self.t_mask])

    def get_all_steps(self):
        return np.copy(self.steps)

    def get_steps_selection(self):
        return np.copy(self.steps[self.t_mask])

    def get_num_steps(self):
        return self.num_steps

    def _get_solution(self):
        if self.p_mask is None or self.t_mask is None or self.variable is None:
            self.solution = None
            return
        self.solution = self.complete_solution[self.t_mask][:, self.p_mask]
        if self.variable == 'coordinates':
            initial_position = self.initial_coordinates[self.p_mask]
            self.solution += initial_position
        if self.component is not None:
            self.solution = self.solution[:, :, self.component]

    def reset_points_selection(self):
        self.select_points(np.full(self.complete_size, True))

    def select_points(self, p_mask):
        # check mask
        if p_mask.shape != (self.complete_size,):
            raise ValueError(f'Point mask has shape {p_mask.shape}, expected shape {(self.complete_size,)}')
        if np.sum(p_mask) == 0:
            tools.print_info('Empty points mask provided', layout='warning')
        self.p_mask = p_mask
        self.size = np.sum(self.p_mask)

    def reset_times_selection(self):
        self.select_times(np.full(self.complete_num_steps + 1, True))

    def select_times(self, t_mask):
        # check mask
        if t_mask.shape != (self.complete_num_steps + 1,):
            raise ValueError(f'Point mask has shape {t_mask.shape}, expected shape {(self.complete_num_steps + 1,)}')
        elif np.sum(t_mask) == 0:
            tools.print_info('Empty time mask provided', layout='warning')
        self.t_mask = t_mask
        self.num_steps = np.sum(self.t_mask) - 1

    def get_values(self, variable=None, component=None):
        self._set_variable(variable, component=component)
        if self.variable is None:
            raise ValueError(f'Specify variable: {", ".join(self.available_vars)}')
        self._get_solution()
        return np.copy(self.solution)

    def __repr__(self):
        return f'SubSet with {self.size} point{"" if self.size == 1 else "s"} and {self.num_steps} ' \
               f'time step{"" if self.num_steps == 1 else "s"} of {self.run}'


class BaseFigure:
    def __init__(self, figure_name=None, aspect='auto'):
        # check figure name, unique to figure
        if figure_name is not None:
            if not plt.fignum_exists(figure_name):
                self.figure_name = figure_name
            else:
                self.figure_name = self._check_figure_name(figure_name)
                tools.print_info(f'Figure_name "{figure_name}" already exists, changed to "{self.figure_name}"',
                                 layout='warning')
        else:
            figure_name = f'{self.__class__.__name__}'
            self.figure_name = self._check_figure_name(figure_name)

        self.subsets = []
        self.aspect = aspect
        self.figure = plt.figure(self.figure_name)
        self.figure_number = self.figure.number
        self.ax = None

        super().__init__()

    @staticmethod
    def _check_figure_name(figure_name):
        i = 1
        while plt.fignum_exists(f'{figure_name} {i}'):
            i += 1
        return f'{figure_name} {i}'

    def _check_subset_name(self, subset, name):
        if name is None:
            name = f'{subset.run.case_name} {subset.model_part_name}'
        if name in self.get_subset_names():
            i = 1
            while f'{name} {i}' in self.get_subset_names():
                i += 1
            tools.print_info(f'SubSet name "{name}" already present, changed to "{name} {i}"')
            name = f'{name} {i}'
        return name

    def add_subset(self, subset, name=None):
        name = self._check_subset_name(subset, name)
        self.subsets.append((dict(name=name, subset=subset)))

    def remove_subset(self, name_or_index):
        self.subsets.pop(self._get_index(name_or_index))
        self.update_figure_layout()

    def get_subsets(self):
        return [subset['subset'] for subset in self.subsets]

    def _get_index(self, name_or_index):
        names = self.get_subset_names()
        if isinstance(name_or_index, int) and name_or_index in range(len(self.subsets)):
            index = name_or_index
        elif isinstance(name_or_index, str) and name_or_index in names:
            index = names.index(name_or_index)
        else:
            raise ValueError(f'Unknown value "{name_or_index}", provide index of SubSet to remove or its name\n'
                             f'Available options:\n\t{", ".join(names)}')
        return index

    def _get_subset_property(self, name_or_index, prop):
        return self.subsets[self._get_index(name_or_index)][prop]

    def get_subset(self, name_or_index):
        return self._get_subset_property(name_or_index, 'subset')

    def get_subset_names(self):
        return [subset['name'] for subset in self.subsets]

    def update_figure_layout(self):
        if self.ax is not None:
            try:  # 'auto' is not implemented for Axes3D for earlier Python versions
                self.ax.set_aspect(self.aspect, adjustable='box')
            except NotImplementedError:
                self.ax.set_box_aspect([1, 1, 1])
            if len(self.subsets) > 1:
                self.ax.legend()
        self.figure.tight_layout()

    def get_figure(self):
        return self.figure

    def get_figure_name(self):
        return self.figure_name

    def get_ax(self):
        return self.ax

    def make_active(self):
        plt.figure(self.figure_name)

    def save(self, file_path=None):
        if file_path is None:
            file_path = f'{self.figure_name}.png'
        self.figure.savefig(file_path)


class Figure(BaseFigure):
    def __init__(self, figure_name=None, aspect='auto', print_function=None, text_box_style=None):
        self.print_text = (print_function is not False)
        if self.print_text:
            self.text = None
            self.print_function = self.print_time if print_function is None else print_function
            self.text_box_style = text_box_style if text_box_style is not None else \
                dict(facecolor='silver', edgecolor='black', pad=5.0, alpha=0.5)

        self.base_time_step = None  # potentially only used for initialization
        self.time = None  # potentially only used for initialization
        self.base_dt = None  # common minimal time step size
        self.base_ts_start = None  # first base time step
        self.base_ts_end = None  # final base time step

        super().__init__(figure_name=figure_name, aspect=aspect)

    def _add_first_subset(self, subset, *args, **kwargs):
        if self.base_dt is None:  # first subset
            self.base_dt = subset.run.dt

        # set time and time step for initialization
        self.base_time_step = subset.get_steps_selection()[0]
        self.time = self.base_time_step * self.base_dt

    def add_subset(self, subset, name=None):
        super().add_subset(subset, name=name)
        self._update_time_step(subset)

    def remove_subset(self, name_or_index):
        super().remove_subset(name_or_index)
        self._reset_time_steps()

    def _update_time_step(self, subset):
        old_base_dt = self.base_dt

        # update base dt
        lcm = self.lcm_of_denominator(self.base_dt, subset.run.dt)
        if lcm > 1:
            self.base_dt = 1 / self.lcm_of_denominator(self.base_dt, subset.run.dt)
        elif lcm == 1.0:
            self.base_dt = subset.run.dt
        else:
            self.base_dt = np.gcd(self.base_dt, subset.run.dt)

        if old_base_dt != self.base_dt:
            tools.print_info(f'{self.figure_name}: SubSets with different time steps detected'
                             f'\n\tUsed time step size is {self.base_dt}')

        # update dt_ratio, base_ts_start and base_ts_end
        for ss in self.subsets:
            ss['dt_ratio'] = dt_ratio = int(ss['subset'].run.dt / self.base_dt)  # ratio of subset dt to base dt
            ss['base_ts_start'] = ss['subset'].get_steps_selection()[0] * dt_ratio
            ss['base_ts_end'] = (ss['subset'].get_steps_selection()[0]
                                 + ss['subset'].get_num_steps() + 1) * dt_ratio - 1

        # update global base time step start and end
        self.base_ts_start = min([ss['base_ts_start'] for ss in self.subsets])
        self.base_ts_end = max([ss['base_ts_end'] for ss in self.subsets])

        # update base time step
        self.base_time_step = self._convert_to_time_step(self.time)

    def _reset_time_steps(self):
        self.base_dt = self.subsets[0]['subset'].run.dt
        for ss in self.subsets:
            self._update_time_step(ss['subset'])

    @staticmethod
    # find common denominator for updating time step sizes
    def lcm_of_denominator(a, b, max_denominator=1e9):
        """
        Finds the least common multiple of denominators of two floats that have been converted to fractions.
        :param a: First float.
        :param b: Second float.
        :param max_denominator: Max denominator allowed for conversion to fraction.
        :return: Least common divider of the denominators.
        """
        fraction_a = Fraction(a).limit_denominator(int(max_denominator))
        fraction_b = Fraction(b).limit_denominator(int(max_denominator))
        multiple = np.lcm(fraction_a.denominator, fraction_b.denominator)
        return multiple

    def _convert_to_time_step(self, time, dt=None, max_denominator=1e9):
        if dt is None:
            dt = self.base_dt
        return round(Fraction(time).limit_denominator(int(max_denominator)) / dt)

    @staticmethod
    def print_time(time):
        if time is None:
            return f"time = - s"
        elif time >= 1e-4 or time == 0:
            return f'time = {time:.4f} s'
        else:
            return f'time = {time * 1e6:.3f} µs'

    def initialize(self):
        pass

    def update(self, base_time_step, update_subsets=False):
        pass

    def update_subsets(self):
        self.update(self.base_time_step, update_subsets=True)

    def update_figure_layout(self):
        super().update_figure_layout()
        if self.print_text:
            self.text.set_text(self.print_time(self.base_time_step * self.base_dt))
        self.figure.tight_layout()


class Figure2d(Figure):
    def __init__(self, subset, x_variable, y_variable, x_component=None, y_component=None, time_step=None, time=None,
                 name=None, text_location=(0.1, 0.1), aspect='auto', **kwargs):
        super().__init__(aspect=aspect, **kwargs)

        self.ax = self.figure.add_subplot()
        self.lines = []

        # add text to figure
        if self.print_text:
            self.text = self.ax.text(*text_location, self.print_time(self.time), transform=self.ax.transAxes,
                                     bbox=self.text_box_style)

        # add subsets if provided
        settings = dict(x_component=x_component, y_component=y_component,
                        time_step=time_step, time=time)
        if isinstance(subset, SubSet):
            self.add_subset(subset, x_variable, y_variable, name=name, **settings)
        elif isinstance(subset, collections.abc.Sequence) and all(isinstance(ss, SubSet) for ss in subset):
            if name is None:
                names = (None,) * len(subset)
            elif isinstance(name, str) or name is None:
                names = (f'{name} {i}' for i in range(len(subset)))
            elif isinstance(name, collections.abc.Sequence) and all(isinstance(n, str) for n in name):
                names = name
                if len(name) != len(subset):
                    raise ValueError('"name" should have the same length as the number of arguments provided')
            else:
                raise ValueError('"name" should be None, string or iterator of strings')
            for subset, name in zip(subset, names):
                if isinstance(subset, SubSet):
                    self.add_subset(subset, x_variable, y_variable, name=name, **settings)
        else:
            raise ValueError('Provided argument(s) should be of type SubSet or iterator of SubSets')

    def _add_first_subset(self, subset, *args, x_component=None, y_component=None, time_step=None, time=None,
                          name=None):
        super()._add_first_subset(subset, time_step=time_step, time=time)

        x_variable, y_variable = args

        # check variables
        available_vars = subset.get_available_variables() + ['initial_coordinates']
        if x_variable not in available_vars:
            raise ValueError(f'x_variable "{x_variable}" is unknown, choose from {", ".join(available_vars)}')
        if y_variable not in available_vars:
            raise ValueError(f'y_variable "{y_variable}" is unknown, choose from {", ".join(available_vars)}')
        self.x_variable = x_variable
        self.y_variable = y_variable

        # check components
        if 'coordinates' not in x_variable and variables_dimensions[x_variable] == 1:
            self.x_component = 0
        elif x_component is None:
            raise ValueError(f'"{x_variable}" has {variables_dimensions[x_variable]} components, '
                             f'specify x_component, choose from 0, 1, 2 or "x", "y", "z"')
        elif x_component in ('x', 'y', 'z'):
            self.x_component = ('x', 'y', 'z').index(x_component)
        elif x_component in (0, 1, 2):
            self.x_component = x_component
        else:
            raise ValueError(f'x_component "{x_component}" invalid, choose from 0, 1, 2 or "x", "y", "z"')
        if 'coordinates' not in y_variable and variables_dimensions[y_variable] == 1:
            self.y_component = 0
        elif y_component is None:
            raise ValueError(f'"{y_variable}" has {variables_dimensions[y_variable]} components, '
                             f'specify y_component, choose from 0, 1, 2 or "x", "y", "z"')
        elif y_component in ('x', 'y', 'z'):
            self.y_component = ('x', 'y', 'z').index(y_component)
        elif y_component in (0, 1, 2):
            self.y_component = y_component
        else:
            raise ValueError(f'y_component "{y_component}" invalid, choose from 0, 1, 2 or "x", "y", "z"')

        x_label = self.x_variable if x_component is None else f'{("x", "y", "z")[self.x_component]} {self.x_variable}'
        y_label = self.y_variable if y_component is None else f'{("x", "y", "z")[self.y_component]} {self.y_variable}'
        self.ax.set_xlabel(x_label.replace('_', ' '))
        self.ax.set_ylabel(y_label.replace('_', ' '))

    def add_subset(self, subset, *args, name=None, **kwargs):
        # check if first subset
        if len(self.subsets) == 0:
            self._add_first_subset(subset, *args, **kwargs)

        # check if provided instance is SubSet
        if not isinstance(subset, SubSet):
            raise ValueError('Provided argument should be of type SubSet')

        # check if variables are present in subset
        for variable in (self.x_variable, self.y_variable):
            if variable not in subset.get_available_variables() + ['initial_coordinates']:
                raise ValueError(f'Cannot add subset: "{variable}" not present in provided SubSet: \n{subset}')

        # check if time step is present in subset
        time_step = self._convert_to_time_step(self.time, dt=subset.run.dt)
        if time_step not in subset.steps:
            tools.print_info(f'{self.figure_name}: time "{self.time}" not present in provided SubSet:\n{subset}',
                             layout='warning')

        # add subset
        name = self._check_subset_name(subset, name)

        def min_max(a):
            return np.min(a), np.max(a)

        if self.x_variable == 'initial_coordinates':
            x = x_tot = subset.get_initial_coordinates_selection()[:, self.x_component]
        else:
            x = x_tot = subset.get_values(self.x_variable, component=self.x_component)
            x = x[time_step - subset.get_steps_selection()[0]]
        if self.y_variable == 'initial_coordinates':
            y = y_tot = subset.get_initial_coordinates_selection()[:, self.y_component]
        else:
            y = y_tot = subset.get_values(self.y_variable, component=self.y_component)
            y = y[time_step - subset.get_steps_selection()[0]]

        x_min, x_max = min_max(x_tot)
        y_min, y_max = min_max(y_tot)

        line, = self.ax.plot(x, y, label=name)

        self.subsets.append(dict(name=name, subset=subset, line=line, x=x_tot, y=y_tot, x_limits=[x_min, x_max],
                                 y_limits=[y_min, y_max]))
        self._update_time_step(subset)

        self.update_figure_layout()

    def remove_subset(self, name_or_index):
        self.get_line(name_or_index).remove()
        super().remove_subset(name_or_index)

    def update_figure_layout(self):
        self.update_limits()
        super().update_figure_layout()

    def update_limits(self):
        self.update_x_limits()
        self.update_y_limits()

    def _update_limits(self, limits_name):
        minimum, maximum = np.inf, -np.inf
        for subset in self.subsets:
            limits = subset[limits_name]
            minimum = min(minimum, limits[0])
            maximum = max(maximum, limits[1])
        return minimum, maximum

    def update_x_limits(self):
        min_value, max_value = self._update_limits('x_limits')
        if min_value < max_value:
            margin = (max_value - min_value) * 0.05
            self.ax.set_xlim([min_value - margin, max_value + margin])

    def update_y_limits(self):
        min_value, max_value = self._update_limits('y_limits')
        if min_value < max_value:
            margin = (max_value - min_value) * 0.05
            self.ax.set_ylim([min_value - margin, max_value + margin])

    def get_lines(self):
        return [subset['line'] for subset in self.subsets]

    def get_line(self, name_or_index):
        return self._get_subset_property(name_or_index, 'line')

    def initialize(self):
        for subset in self.subsets:
            self.empty_line(subset['line'])

    @staticmethod
    def empty_line(line):
        nan_array = np.empty(line.get_data()[0].shape)
        nan_array.fill(np.nan)
        line.set_ydata(nan_array)

    def update(self, base_time_step, update_subsets=False):
        for ss in self.subsets:
            subset, line = ss['subset'], ss['line']
            base_ts_start, base_ts_end, dt_ratio = ss['base_ts_start'], ss['base_ts_end'], ss['dt_ratio']
            if not (base_ts_start <= base_time_step <= base_ts_end):
                self.empty_line(line)
            else:
                time_index = base_time_step // dt_ratio - subset.get_steps_selection()[0]
                # update x-variable
                if self.x_variable != 'initial_coordinates':
                    if update_subsets:
                        ss['x'] = subset.get_values(self.x_variable, component=self.x_component)
                    x = ss['x'][time_index]
                    line.set_xdata(x)

                # update y-variable
                if self.y_variable != 'initial_coordinates':
                    if update_subsets:
                        ss['y'] = subset.get_values(self.y_variable, component=self.y_component)
                    y = ss['y'][time_index]
                    line.set_ydata(y)

        if self.print_text:
            self.text.set_text(self.print_time(base_time_step * self.base_dt))

    def __repr__(self):
        x = self.x_variable if self.x_component is None else f'{("x", "y", "z")[self.x_component]} {self.x_variable}'
        y = self.y_variable if self.y_component is None else f'{("x", "y", "z")[self.y_component]} {self.y_variable}'
        return f'{self.figure_name} named {self.figure_name} showing ' \
               f'{y.replace("_", " ")} versus {x.replace("_", " ")} of\n{[ss["subset"] for ss in self.subsets]}'


class Figure3d(Figure):
    def __init__(self, subset, variable=None, color_by=None, component=None, time_step=None, time=None, name=None,
                 text_location=(0, 0), aspect='equal', **kwargs):
        super().__init__(aspect=aspect, **kwargs)

        self.variable = None
        self.color_by = None
        self.component = None

        plt.figure(self.figure_number)  # make figure active
        self.ax = plt.axes(projection='3d')
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')
        self.color_bar = None

        # add text to figure
        if self.print_text:
            self.text = self.ax.text2D(*text_location, self.print_time(self.time), transform=self.ax.transAxes,
                                       bbox=self.text_box_style)

        # add subsets if provided
        settings = dict(variable=variable, color_by=color_by, component=component, time_step=time_step, time=time)
        if isinstance(subset, SubSet):
            self.add_subset(subset, name=name, **settings)
        elif isinstance(subset, collections.abc.Sequence) and all(isinstance(ss, SubSet) for ss in subset):
            if name is None:
                names = (None,) * len(subset)
            elif isinstance(name, str) or name is None:
                names = (f'{name} {i}' for i in range(len(subset)))
            elif isinstance(name, collections.abc.Sequence) and all(isinstance(n, str) for n in name):
                names = name
                if len(name) != len(subset):
                    raise ValueError('"name" should have the same length as the number of arguments provided')
            else:
                raise ValueError('"name" should be None, string or iterator of strings')
            for subset, name in zip(subset, names):
                if isinstance(subset, SubSet):
                    self.add_subset(subset, name=name, **settings)
        else:
            raise ValueError('Provided argument(s) should be of type SubSet or iterator of SubSets')

    def _add_first_subset(self, subset, variable=None, color_by=None, component=None, time_step=None, time=None):
        super()._add_first_subset(subset, time_step=time_step, time=time)

        # check variable
        if variable is None:
            self.variable = 'coordinates' if 'displacement' in subset.available_vars else 'initial_coordinates'
        elif variable in ('coordinates', 'initial_coordinates'):
            self.variable = variable
        else:
            available_vars = subset.get_available_variables() + ['initial_coordinates']
            raise ValueError(f'Variable "{variable}"s unknown, choose from {", ".join(available_vars)}')

        # check color_by
        if self.color_by is None:
            # check color_by
            if color_by is not None and color_by not in subset.get_available_variables():
                raise ValueError(
                    f'color_by "{color_by}" is unknown, choose from {", ".join(subset.get_available_variables())}')
            else:
                self.color_by = color_by

        # check component
        if self.component is None:
            # check component
            if component in (0, 1, 2, None):
                self.component = component
            elif component in ('x', 'y', 'z'):
                self.component = ('x', 'y', 'z').index(component)
            else:
                raise ValueError(f'component "{component}" invalid, choose from 0, 1, 2 or "x", "y", "z" or None')

    def add_subset(self, subset, name=None, cmap=plt.cm.coolwarm, **kwargs):
        # check if first subset
        if len(self.subsets) == 0:
            self._add_first_subset(subset, **kwargs)

        # check if provided instance is SubSet
        if not isinstance(subset, SubSet):
            raise ValueError('Provided argument should be of type SubSet')

        # check if variable is present in subset
        if self.variable not in subset.get_available_variables() + ['initial_coordinates']:
            raise ValueError(f'Cannot add subset: "{self.variable}" not present in provided SubSet: \n{subset}')

        # check color_by
        if self.color_by is not None and self.color_by not in subset.get_available_variables():
            raise ValueError(f'Cannot add subset: "{self.color_by}" is not present in provided SubSet:\n{subset}')

        # check if time step is present in subset
        time_step = self._convert_to_time_step(self.time, dt=subset.run.dt)
        if time_step not in subset.steps:
            tools.print_info(f'{self.figure_name}: time "{self.time}" not present in provided SubSet:\n{subset}',
                             layout='warning')

        # add subset
        name = self._check_subset_name(subset, name)

        if self.color_by is None:
            color = color_tot = None
            color_min, color_max = 0, 0
        else:
            color = subset.get_values(self.color_by, component=self.component)
            if len(color.shape) > 2:
                color = np.linalg.norm(color, axis=-1)
            color_min, color_max = np.min(color), np.max(color)
            color_tot = color
            color = color[time_step - subset.get_steps_selection()[0]]

        if self.variable == 'coordinates' or (
                self.variable is None and 'coordinates' in subset.get_available_variables()):
            coordinates = coord_tot = subset.get_values(self.variable)
            coordinates = coordinates[time_step - subset.get_steps_selection()[0]].T
            tools.print_info(f'{self.figure_name}: Instantaneous coordinates used')
        else:
            coordinates = subset.get_initial_coordinates_selection().T
            coord_tot = None
            tools.print_info(f'{self.figure_name}: Initial coordinates used')

        if color is None:
            scatter = self.ax.scatter(*coordinates, label=name)
        elif self.color_bar is None:
            scatter = self.ax.scatter(*coordinates, c=color, label=name, cmap=cmap)
        else:  # make sure all subsets are plotted using the same color scale
            scatter = self.ax.scatter(*coordinates, c=color, label=name, cmap=cmap, norm=self.color_bar.norm)

        self.subsets.append(dict(name=name, subset=subset, scatter=scatter, color_limits=[color_min, color_max],
                            color=color_tot, coordinates=coord_tot))
        self._update_time_step(subset)

        if self.color_bar is None and self.color_by is not None:
            self.color_bar = self.figure.colorbar(scatter, format="%4.2g")
            self.color_bar.set_label(self._repr_component() + self.color_by)

        self.update_figure_layout()

    def remove_subset(self, name_or_index):
        self.get_scatter(name_or_index).remove()
        super().remove_subset(name_or_index)

    def update_figure_layout(self):
        if self.color_by is not None:
            self.update_color_limits()
            self.color_bar.set_label(self._repr_component() + self.color_by)
        super().update_figure_layout()

    def update_color_limits(self):
        minimum, maximum = np.inf, -np.inf
        for subset in self.subsets:
            limits = subset['color_limits']
            minimum = min(minimum, limits[0])
            maximum = max(maximum, limits[1])
        self.color_bar.norm.autoscale((minimum, maximum))

    def get_scatters(self):
        return [subset['scatter'] for subset in self.subsets]

    def get_scatter(self, name_or_index):
        return self._get_subset_property(name_or_index, 'scatter')

    def initialize(self):
        for subset in self.subsets:
            self.empty_scatter(subset['scatter'])

    def empty_scatter(self, scatter):
        if not self.variable == 'initial_coordinates':
            nan_array = np.empty(scatter._offsets3d[0].shape)
            nan_array.fill(np.nan)
            scatter._offsets3d = 3 * (nan_array,)

    def update(self, base_time_step, update_subsets=False):
        for ss in self.subsets:
            subset, scatter = ss['subset'], ss['scatter']
            base_ts_start, base_ts_end, dt_ratio = ss['base_ts_start'], ss['base_ts_end'], ss['dt_ratio']
            if not (base_ts_start <= base_time_step <= base_ts_end):
                self.empty_scatter(scatter)
            else:
                time_index = base_time_step // dt_ratio - subset.get_steps_selection()[0]
                # update coordinates
                if self.variable == 'coordinates' or (
                        self.variable is None and 'coordinates' in subset.get_available_variables()):
                    if update_subsets:
                        ss['coordinates'] = subset.get_values('coordinates')
                    coordinates = ss['coordinates'][time_index].T
                    scatter._offsets3d = tuple(coordinates)
                # update color
                if self.color_by is not None:
                    if update_subsets:
                        color = subset.get_values(self.color_by, component=self.component)
                        if len(color.shape) > 2:
                            color = np.linalg.norm(color, axis=-1)
                        ss['color'] = color
                    color = ss['color'][time_index]

                    scatter.set_array(color)

        if self.print_text:
            self.text.set_text(self.print_time(base_time_step * self.base_dt))

    def _repr_component(self):
        if self.variable is not None and self.color_by is not None and \
                (self.color_by == 'coordinates' or variables_dimensions[self.color_by] == 3):
            if self.component is not None:
                return f'{("x", "y", "z")[self.component]}-component of '
            else:
                return 'magnitude of '
        else:
            return ''

    def __repr__(self):
        return f'{self.figure_name} named {self.figure_name} showing ' \
               f'{self.variable.replace("_", " ")} of \n{[ss["subset"] for ss in self.subsets]}' \
               f'\ncolored by {self._repr_component()}{self.color_by}'


class Plot(Figure):
    def _add_first_subset(self, subset, time_step=None, time=None, **kwargs):
        super()._add_first_subset(subset, **kwargs)

        # check time step input
        if time_step is not None:
            if time is not None:
                tools.print_info(f'{self.figure_name}: both "time_step" and "time" provided, "time" ignored',
                                 layout='warning')
            self.base_time_step = time_step
            self.time = time_step * self.base_dt
        elif time is not None:
            self.base_time_step = self._convert_to_time_step(time)
            self.time = time
        elif time_step is None and subset.get_num_steps() == 0:
            self.base_time_step = subset.get_steps_selection()[0]
            self.time = self.base_time_step * self.base_dt
        else:
            tools.print_info(f'{self.figure_name}: no time or time step provided, '
                             f'initial time step used', layout='warning')

    def set_time(self, time):
        self.time = time
        self.base_time_step = self._convert_to_time_step(self.time)
        self.update(self.base_time_step)
        self.update_figure_layout()

    def set_time_step(self, base_time_step):
        self.base_time_step = base_time_step
        self.time = self.base_time_step * self.base_dt
        tools.print_info(f'{self.figure_name}: Provided time step corresponds to time {self.time}s')
        self.update(self.base_time_step)
        self.update_figure_layout()

    def __repr__(self):
        return super().__repr__() + f'add time {self.print_time(self.base_time_step * self.base_dt)}'


class Plot2d(Figure2d, Plot):
    pass


class Plot3d(Figure3d, Plot):
    pass


class Animation(Figure):
    def __init__(self, *args, func_animation_settings=None, **kwargs):
        self.animation = None

        self.writer = ani.PillowWriter(fps=15, metadata=dict(artist='CoCoNuT'), bitrate=1800)

        self.frames = None
        self.skip = 0
        self.save_count = None
        self.animation = None

        if func_animation_settings is None:
            func_animation_settings = {}
        func_animation_settings = self.set_func_animation_settings(**func_animation_settings)

        super().__init__(*args, **kwargs)

        self.make_animation(**func_animation_settings)

    def set_func_animation_settings(self, **func_animation_settings):
        func_animation_settings = dict(func_animation_settings)  # creates copy
        # frames: iterable, int, generator function, or None, optional
        #     (None) Plots all non-empty frames.
        #     (int) Number of frames (<= number of time steps + 1).
        #     (iterable) Frames to plot.
        # skip: int, default: 0
        #     Number of frames skipped between each shown frame (used if frames is None or int).
        # save_count: int, default: number of frames
        #     Number of frames to cache.
        if 'frames' in func_animation_settings:
            self.frames = func_animation_settings.pop('frames')
        if 'skip' in func_animation_settings:
            self.skip = func_animation_settings.pop('skip')
        if 'save_count' in func_animation_settings:
            self.save_count = func_animation_settings.pop('save_count')

        if self.animation is not None:  # overwrite animation
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                del self.animation
            self.make_animation(**func_animation_settings)
        return func_animation_settings

    def _update_time_step(self, subset):
        old_base_dt = self.base_dt

        super()._update_time_step(subset)

        dt_ratio = int(old_base_dt / self.base_dt)

        self.skip *= dt_ratio
        if isinstance(self.frames, int):
            self.frames *= dt_ratio
        elif isinstance(self.frames, range):
            self.frames = range(self.frames.start * dt_ratio, self.frames.stop * dt_ratio, self.frames.step * dt_ratio)

    def _reset_time_steps(self):
        dt_ratio = self.base_dt / self.subsets[0]['subset'].run.dt

        self.skip = int(self.skip * dt_ratio)
        if isinstance(self.frames, int):
            self.frames = int(self.frames * dt_ratio)
        elif isinstance(self.frames, range):
            self.frames = range(int(self.frames.start * dt_ratio), int(self.frames.stop * dt_ratio),
                                int(self.frames.step * dt_ratio))

        super()._reset_time_steps()

    def make_animation(self, interval=200, repeat=True, blit=False, **kwargs):
        """
        Creates the animation by repeatedly calling the update() function.

        :param interval: int, default: 200
            Delay between frames in milliseconds.
        :param repeat: bool, default: True
            Whether the animation repeats when the sequence of frames is completed.
        :param blit: bool, default: False
            Whether blitting is used to optimize drawing.
        :param kwargs: optional
            Other keyword arguments to passed to the class FuncAnimation()
        :return: FuncAnimation instance
        """
        self.animation = ani.FuncAnimation(self.figure, self.update, frames=self.determine_frames,
                                           save_count=self.determine_save_count(), init_func=self.initialize,
                                           interval=interval, repeat=repeat, blit=blit, **kwargs)

    def determine_save_count(self):
        num_frames_avail = self.base_ts_end - self.base_ts_start + 1
        if self.frames is None:
            return num_frames_avail
        elif isinstance(self.frames, int):
            return self.frames
        else:
            return len(self.frames)

    def determine_frames(self):
        num_frames_avail = self.base_ts_end - self.base_ts_start + 1
        if self.frames is None:
            if self.skip >= num_frames_avail:
                raise ValueError(f'Invalid parameter "skip": {self.skip}\n'
                                 f'value should be lower than maximum number of frames is {num_frames_avail} '
                                 f'(number of base time steps + 1)')
            frames = range(self.base_ts_start, self.base_ts_end + 1, self.skip + 1)
        elif isinstance(self.frames, int):
            if not 0 < self.frames <= num_frames_avail:
                raise ValueError(f'Invalid parameter "frames": {self.frames}\n'
                                 f'maximum number of frames is {num_frames_avail} (number of base time steps + 1),\n'
                                 f'with time step size {self.base_dt} and starting time step {self.base_ts_start}')
            if self.skip >= num_frames_avail:
                raise ValueError(f'Invalid parameter "skip": {self.skip}\n'
                                 f'value should be lower than maximum number of frames is {num_frames_avail} '
                                 f'(number of base time steps + 1)')
            frames = range(self.base_ts_start, self.base_ts_start + self.frames, self.skip + 1)
        elif min(self.frames) < self.base_ts_start or max(self.frames) > self.base_ts_end:
            raise Exception(f'Invalid parameter "frames":\n'
                            f'minimum time step is {self.base_ts_start}, maximum time step is {self.base_ts_end},\n'
                            f"with time step size {self.base_dt} and starting time step {self.base_ts_start}")
        else:
            frames = self.frames
        return iter(frames)

    def get_animation(self):
        return self.animation

    def save(self, file_path=None):
        if file_path is None:
            file_path = f'{self.figure_name}.gif'
        self.animation.save(file_path, writer=self.writer)

    def set_writer(self, writer):
        if not isinstance(writer, ani.AbstractMovieWriter):
            raise ValueError('Provided writer is no AbstractMovieWriter')
        self.writer = writer


class Animation2d(Animation, Figure2d):
    pass


class Animation3d(Animation, Figure3d):
    pass


class TimeEvolution(BaseFigure):
    def __init__(self, subset, variable, component=None, name=None, **kwargs):
        super().__init__(**kwargs)

        self.variable = None

        self.ax = self.figure.add_subplot()
        self.lines = []

        # add subsets if provided
        settings = dict(component=component)
        if isinstance(subset, SubSet):
            self.add_subset(subset, variable, name=name, **settings)
        elif isinstance(subset, collections.abc.Sequence) and all(isinstance(ss, SubSet) for ss in subset):
            if name is None:
                names = (None,) * len(subset)
            elif isinstance(name, str) or name is None:
                names = (f'{name} {i}' for i in range(len(subset)))
            elif isinstance(name, collections.abc.Sequence) and all(isinstance(n, str) for n in name):
                names = name
                if len(name) != len(subset):
                    raise ValueError('"name" should have the same length as the number of arguments provided')
            else:
                raise ValueError('"name" should be None, string or iterator of strings')
            for subset, name in zip(subset, names):
                if isinstance(subset, SubSet):
                    self.add_subset(subset, variable, name=name, **settings)
        else:
            raise ValueError('Provided argument(s) should be of type SubSet or iterator of SubSets')

    def _add_first_subset(self, subset, variable):
        # check variable
        if variable not in subset.get_available_variables():
            raise ValueError(f'x_variable "{variable}" is unknown, choose from '
                             f'{", ".join(subset.get_available_variables())}')
        self.variable = variable

        self.ax.set_xlabel('time')
        self.ax.set_ylabel(self.variable.replace('_', ' '))

    def add_subset(self, subset, *args, component=None, name=None, **kwargs):
        # check if first subset
        if len(self.subsets) == 0:
            self._add_first_subset(subset, *args)

        # check if provided instance is SubSet
        if not isinstance(subset, SubSet):
            raise ValueError('Provided argument should be of type SubSet')

        # check if variables are present in subset
        if self.variable not in subset.get_available_variables():
            raise ValueError(f'Cannot add subset: "{self.variable}" not present in provided SubSet: \n{subset}')

        # check number of points in subset
        if subset.size > 1:
            raise ValueError(f'More than one point in SubSet')

        # check component
        if component not in (0, 1, 2, None):
            if component in ('x', 'y', 'z'):
                component = ('x', 'y', 'z').index(component)
            else:
                raise ValueError(f'component "{component}" invalid, choose from 0, 1, 2 or "x", "y", "z" or None')

        # add subset
        name = self._check_subset_name(subset, name)
        if self.variable == 'coordinates' or variables_dimensions[self.variable] == 3:
            if component is not None:
                name += f' {("x", "y", "z")[component]}-component'
            else:
                name += ' magnitude '

        def min_max(a):
            return np.min(a), np.max(a)

        x = subset.get_times_selection()
        x_min, x_max = min_max(x)

        y = subset.get_values(self.variable, component=component)

        if len(y.shape) == 3:  # calculate magnitude
            y = np.linalg.norm(y, axis=-1)
        y_min, y_max = min_max(y)

        line, = self.ax.plot(x, y, label=name)

        self.subsets.append(dict(name=name, subset=subset, line=line, x_limits=[x_min, x_max], y_limits=[y_min, y_max],
                                 component=component))

        self.update_figure_layout()

    def _repr_component(self, component):
        if self.variable is not None and \
                (self.variable == 'coordinates' or variables_dimensions[self.variable] == 3):
            if component is not None:
                return f'{("x", "y", "z")[component]}-component of '
            else:
                return 'magnitude of '
        else:
            return ''

    def __repr__(self):
        ss = [f'{self._repr_component(ss["component"])}{self.variable.replace("_", " ")} '
              f'in function of time for \n{ss["subset"]}' for ss in self.subsets]
        return f'{self.figure_name} named {self.figure_name} showing:\n' \
               f'{", ".join(ss)}'


def accept_post_process(subset, interface_name, kwargs):
    if isinstance(subset, PostProcess):
        post_process = subset
        subset = [post_process.add_subset(interface=interface_name, model_part=mp) for
                  mp in post_process.get_model_part_names(interface_name)]
        kwargs.update({'name': post_process.get_model_part_names(interface_name)})
    return subset, kwargs


class Plot2dDisplacement(Plot2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, 'initial_coordinates', 'displacement', x_component=x_component,
                         y_component=y_component, **kwargs)


class Plot2dCoordinates(Plot2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, 'coordinates', 'coordinates', x_component=x_component, y_component=y_component,
                         aspect='equal', **kwargs)


class Plot2dPressure(Plot2d):
    def __init__(self, subset, x_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, 'initial_coordinates', 'pressure', x_component=x_component, **kwargs)


class Plot2dTraction(Plot2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, 'initial_coordinates', 'traction', x_component=x_component,
                         y_component=y_component, **kwargs)


class Plot3dDisplacement(Plot3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, color_by='displacement', **kwargs)


class Plot3dCoordinates(Plot3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, **kwargs)


class Plot3dPressure(Plot3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, color_by='pressure', **kwargs)


class Plot3dTraction(Plot3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, color_by='traction', **kwargs)


class Animation2dDisplacement(Animation2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, 'initial_coordinates', 'displacement', x_component=x_component,
                         y_component=y_component, **kwargs)


class Animation2dCoordinates(Animation2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, 'coordinates', 'coordinates', x_component=x_component, y_component=y_component,
                         aspect='equal', **kwargs)


class Animation2dPressure(Animation2d):
    def __init__(self, subset, x_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, 'initial_coordinates', 'pressure', x_component=x_component, **kwargs)


class Animation2dTraction(Animation2d):
    def __init__(self, subset, x_component, y_component, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, 'initial_coordinates', 'traction', x_component=x_component,
                         y_component=y_component, **kwargs)


class Animation3dDisplacement(Animation3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, color_by='displacement', **kwargs)


class Animation3dCoordinates(Animation3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_x', kwargs)
        super().__init__(subset, **kwargs)


class Animation3dPressure(Animation3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, color_by='pressure', **kwargs)


class Animation3dTraction(Animation3d):
    def __init__(self, subset, **kwargs):
        subset, kwargs = accept_post_process(subset, 'interface_y', kwargs)
        super().__init__(subset, color_by='traction', **kwargs)
