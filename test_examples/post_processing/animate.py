from coconut import data_structure
from coconut.coupling_components import tools

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pickle

# This files contains the calss Animate, which can be used make an animation of different cases or a plot at a certain
# time step.
# To generate a result file, include a boolean "save_results" in the settings of the coupled solver with value True.
# Give a name to the case by including the string {"name": "a_name"} in the settings of the coupled solver .


class Animate:
    animations_list = []
    min = None
    max = None

    def __init__(self, figure, solution, interface, model_part=None, variable=None, name=None):
        """
        Creates Animate instance
            self.info: list e.g. [("mp_a", ["PRESSURE", "TRACTION"]), ("mp_b, "DISPLACEMENT")]
                this dictates the order of the values in solution and coordinates:
                e.g. p0, p1, ..., pm, tx0, ty0, tz0, tx1, ty1, tz1,, ..., txm, tym, tzm, dx0, dy0, dz0, ...,
                    dxm, dym, dzm
                where p and t are pressure and traction on mp_a and d is displacement on mp_b
            self.coordinates: np.array contains the initial coordinates of the nodes on the interface
                as given in self.info
                  e.g. x0, y0, z0, x1, y1, z1, ...

        :param figure: (plt.figure) figure where animation is created
        :param solution : (np.array) contains as columns the solution of each time step, order is dictated by self.info
        :param interface: (Interface) object interface to be plotted
        :param model_part: (string) model part to be plotted (optional if there is only one)
        :param variable: (string) variable to be plotted (optional if there is only one corresponding to the model part)
        :param name: (string) the name in the legend (optional)
        """
        """"

        """
        self.fig = figure
        self.complete_solution = solution
        self.interface = interface
        self.info = interface.model_parts_variables
        self.coordinates = interface.GetInitialCoordinates()
        self.m = int(self.coordinates.size / 3)  # number of nodes
        self.time_steps = solution.shape[1] - 1
        self.animation = None
        self.mask = None
        self.argsort = None
        self.absicissa = None
        self.solution = None
        self.displacement = None
        self.line = None

        # check that model_part or variables are given if not unique
        if model_part is None:
            if len(self.info) > 1:
                raise Exception(f"Specify model_part: more than one present: {self.info}")
            else:
                model_part_name = self.info[0][0]
                index = 0
        elif model_part not in [t[0] for t in self.info]:
            raise Exception(f"Given model_part '{model_part}' is not found")
        else:
            model_part_name = model_part
            index = self.info[self.info.find(model_part)]
        if variable is None:
            if len(self.info[index][1].list()) > 1:
                raise Exception(f"Specify variable: more than one present: {self.info[index][1].list()}")
            else:
                variable_name = self.info[index][1].list()[0].GetString()
        else:
            variable = variable.upper()
            if variable not in [var.GetString() for var in self.info[index][1].list()]:
                raise Exception(f"Given variable '{variable}' is not found")
            else:
                variable_name = variable

        self.name = name if name is not None else model_part_name + ': ' + variable_name

        # find location of data
        index = 0
        self.displacement_available = False
        for mp_name, var_names in self.info:
            mp = interface.model.GetModelPart(mp_name)
            number_of_nodes = sum(1 for _ in mp.Nodes)
            for var_name in [var.GetString() for var in var_names.list()]:
                correct_location = True if var_name == variable_name and mp_name == model_part_name else False
                displacement_location = True if var_name == "DISPLACEMENT" and mp_name == model_part_name else False
                var = vars(data_structure)[var_name]
                if var.Type() is "Double":
                    dimension = 1
                elif var.Type() is "Array":
                    dimension = 3
                else:
                    raise NotImplementedError('Only "Double" and "Array" Variables implemented.')
                if correct_location:
                    self.start_index = index
                    self.dimension = dimension
                    self.end_index = index + self.dimension * self.m
                    if self.m != number_of_nodes:
                        raise Exception("Number of coordinates do not match.")
                if displacement_location:
                    self.start_index_displacement = index
                    self.end_index_displacement = index + 3 * self.m
                    self.displacement_available = True
                index += dimension * number_of_nodes

        if not self.displacement_available:
            tools.Print(f"{self.name}: Nodes positions are not updated, because no 'DISPLACEMENT' available.",
                        layout='warning')

        if index != self.complete_solution.shape[0]:
            raise Exception("Size of provided solution data does not match interface.")

        # add animation instance to class list
        self.animations_list.append(self)

    def Initialize(self, mask_x, mask_y, mask_z, absicissa, component):
        """
        This method selects which points to plot and how to sort them.

        :param mask_x: (ndarray) selects points based on x-coordinate
        :param mask_y: (ndarray) selects points based on y-coordinate
        :param mask_z: (ndarray) selects points based on z-coordinate
        :param absicissa: (int) absicissa direction: 0 for x-axis, 1 for y-axis and 2 for z-axis
        :param component: (int) which component to plot if variable is vector:
                                0 for x-axis, 1 for y-axis and 2 for z-axis
        """
        # chose which nodes to plot
        self.mask = mask_x & mask_y & mask_z

        # chose sort direction
        self.argsort = np.argsort(self.coordinates[absicissa::3][self.mask])
        self.absicissa = self.coordinates[absicissa::3][self.mask][self.argsort]

        if self.dimension == 1:
            component = 0
        self.solution = [self.select(self.complete_solution[:, i], component) for i in range(self.time_steps + 1)]
        if self.displacement_available:
            self.displacement = [self.select_displacement(self.complete_solution[:, i], absicissa)
                                 for i in range(self.time_steps + 1)]

        # adjust scale
        minimum = min([s.min() for s in self.solution])
        maximum = max([s.max() for s in self.solution])
        Animate.min = minimum if Animate.min is None else min(Animate.min, minimum)
        Animate.max = maximum if Animate.max is None else max(Animate.max, maximum)
        margin = (Animate.max - Animate.min) * 0.05
        plt.ylim([Animate.min - margin, Animate.max + margin])

        self.line, = plt.plot(self.absicissa, self.solution[0], label=self.name)

    def select(self, array, component):
        # select correct model_part and variable data and order them
        return array[self.start_index + component: self.end_index: self.dimension][self.mask][self.argsort]

    def select_displacement(self, array, absicissa):
        # select correct model_part and variable node displacement data and order them
        return array[self.start_index_displacement + absicissa: self.end_index_displacement: 3][self.mask][self.argsort]

    def case_init(self):
        self.line.set_ydata([np.nan] * self.absicissa.size)
        return self.line,

    def case_animate(self, ts):
        if self.displacement_available:
            self.line.set_xdata(self.absicissa + self.displacement[ts])
        else:
            self.line.set_xdata(self.absicissa)
        self.line.set_ydata(self.solution[ts])
        return self.line,

    def init(self):  # only required for blitting to give a clean slate.
        lines = ()
        for animation in Animate.animations_list:
            lines += animation.case_init()
        return lines,

    def animate(self, ts):
        lines = ()
        for animation in Animate.animations_list:
            lines += animation.case_animate(ts)
        return lines,

    def MakeAnimation(self, interval=100, blit=False, save_count=100, repeat=True, frames=None):
        # inteval: interval between frames in ms
        # frames: (int) number of frames (<= number of time steps + 1)
        #         (iterable) frames to plot (index <= number of time steps)
        if self.line is None:
            raise Exception("Animate object has not yet been initialized.")
        if frames is None:
            frames = self.time_steps + 1
        elif type(frames) is int:
            if frames > self.time_steps + 1:
                raise Exception(f"Time step out of range: maximum number of frames is {self.time_steps + 1} (number of "
                                f"time steps + 1).")
            else:
                pass
        elif max(frames) > self.time_steps:
            raise Exception(f"Time step out of range: maximum time step is {self.time_steps}.")
        self.animation = ani.FuncAnimation(self.fig, self.animate, init_func=self.init, interval=interval,
                                           blit=blit, save_count=save_count, repeat=repeat, frames=frames)

    def MakePlot(self, time_step):
        # time_step: time step at which plot is made
        if self.line is None:
            raise Exception("Animate object has not yet been initialized.")
        if time_step > self.time_steps:
            raise Exception(f"Time step out of range: maximum time step is {self.time_steps}")
        self.animate(time_step)


# different cases to be plotted
common_path = "../../test_examples/"
case_names = ["results"]
case_paths = ["tube_tube_flow_tube_structure/results"]

# load cases
results = {}
for name, path in zip(case_names, case_paths):
    results.update({name: pickle.load(open(os.path.join(common_path, path), 'rb'))})

# make figure and create animation for each case
figure = plt.figure()
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
line_styles = ['-', '--', ':', '-.']
for j, name in enumerate(case_names):
    solution = results[name]["solution_x"]
    interface = results[name]["interface_x"]
    # create animate object
    animation = Animate(figure, solution, interface, variable="displacement", name=name)
    # select points and component of variable to plot
    coordinates = animation.coordinates

    # example 1: python solver (YZ-plane)
    python_solver = True
    if python_solver:
        mask_x = (coordinates[::3] > -np.inf)
        mask_y = (abs(coordinates[1::3]) > 0)
        mask_z = (coordinates[2::3] > -np.inf)
        absicissa = 2  # z-axis
        component = 1  # y-component

    # example 2: fluent solver (XY-plane)
    fluent = False
    if fluent:
        mask_x = (coordinates[::3] > -np.inf)
        mask_y = (coordinates[1::3] > 0)
        mask_z = (abs(coordinates[2::3]) < 1e-16)
        # mask_z = (coordinates[2::3] > 0) & (coordinates[2::3] < 0.0005)
        absicissa = 0  # x-axis
        component = 1  # y-component

    animation.Initialize(mask_x, mask_y, mask_z, absicissa, component)
    animation.line.set_color(colors[j % len(colors)])
    animation.line.set_linestyle(line_styles[0])
    animation.line.set_marker('o')
    animation.line.set_markersize(2)

plt.ylabel("displacement (m)")
plt.xlabel("axial coordinate (m)")
plt.legend()

animation.MakeAnimation()
# animation.MakePlot(50)

save = False
movie_name = "displacement.mp4"
if save:
    # set up formatting for the movie files: mp4-file
    plt.rcParams['animation.ffmpeg_path'] = u'/apps/SL6.3/FFmpeg/3.1.3/bin/ffmpeg'  # path to ffmpeg conversion tool
    writer = ani.FFMpegFileWriter(codec='mpeg1video', metadata=dict(artist='NicolasDelaissé'), fps=24, bitrate=2000)

    animation.animation.save(movie_name, writer=writer)

plt.show()
