from coconut import solver_modules

import time
from contextlib import contextmanager
import numpy as np
import warnings
import os
from os.path import join
import subprocess
import pickle
import importlib.util
import shutil


def create_instance(settings, if_not_defined=None):
    if 'type' not in settings and if_not_defined is not None:
        # create dummy object denoted by if_not_defined
        object_type = if_not_defined
    else:
        # create instance of given class based on settings dict
        object_type = settings['type']
    object_module = __import__('coconut.coupling_components.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


class CocoMessages:

    def __init__(self, working_directory, max_wait_time=1e4, timed_out_action=None, poll_time=0.1):
        self.working_directory = working_directory
        self.max_wait_time = max_wait_time  # in seconds
        self.timed_out_action = timed_out_action  # receives message as argument
        self.poll_time = poll_time  # in seconds
        self.process = None  # process that will be polled

    def set_process(self, process):
        self.process = process

    def send_message(self, message):
        file = join(self.working_directory, message + '.coco')
        open(file, 'w').close()
        return

    def wait_message(self, message):
        cumul_time = 0
        polled = 0
        file = join(self.working_directory, message + '.coco')
        while not os.path.isfile(file):
            time.sleep(0.001)
            cumul_time += 0.001
            if self.process is not None and (cumul_time // self.poll_time) > polled:
                polled = cumul_time // self.poll_time
                if self.process.poll() is not None:  # process has terminated
                    raise RuntimeError(f'Solver process "{self.process.args[:35]}..." '
                                       f'has terminated unexpectedly while waiting for message: {message}.coco')
            elif cumul_time > self.max_wait_time:
                if self.timed_out_action is not None:
                    self.timed_out_action(message)
                raise RuntimeError(f'CoCoNuT timed out, waiting for message: {message}.coco')
        os.remove(file)
        return

    def check_message(self, message):
        file = join(self.working_directory, message + '.coco')
        if os.path.isfile(file):
            os.remove(file)
            return True
        return False

    def remove_all_messages(self):
        for file_name in os.listdir(self.working_directory):
            if file_name.endswith('.coco'):
                file = join(self.working_directory, file_name)
                os.remove(file)


class LayoutStyles:
    styles = {'plain': '\033[0m',
              'bold': '\033[1m',
              'underline': '\033[4m',
              'inverse': '\033[7m',
              'negative': '\033[97m',
              'info': '\033[94m',
              'warning': '\033[91m',
              'fail': '\033[41m',
              'grey': '\033[90m',
              'red': '\033[91m',
              'green': '\033[92m',
              'yellow': '\033[93m',
              'blue': '\033[94m',
              'magenta': '\033[95m',
              'cyan': '\033[96m',
              'white': '\033[1;8m',
              'black': '\033[1;90m',
              }

    def get(self, layout):
        if layout is None:
            return self.styles['plain']
        self.check(layout)
        return self.styles[layout]

    def check(self, layout):
        if layout not in self.styles:
            raise ValueError(f'Layout style "{layout}" is not implemented, correct layout styles are: '
                             'plain, bold, underline, inverse, negative, info, warning, fail, grey, red, green, yellow,'
                             'blue, magenta, cyan, white and black.')


layout_style = LayoutStyles()


# print_info: printing with color
#  @param args          The arguments to be printed
#  @param layout        The layout to be used: plain, bold, underline, inverse, negative, info, warning, fail, grey,
#                       red, green, yellow, blue, magenta, cyan, white, black
def print_info(*args, layout=None, **kwargs):
    if layout is None:
        print(''.join(map(str, args)), **kwargs)
    else:
        print(layout_style.get(layout), ''.join(map(str, args)), layout_style.get('plain'), **kwargs)


# updatePre: update preceding text, used in structure printing
#  @param pre         Preceding text ending with '├─' or '└─'
def update_pre(pre):
    """
    This function is used to update the preceding string in a structured tree-like format as such:
    Name
    ├─Name
    │ XXX
    └─Name
      XXX
    """
    if pre[-2:] == '├─':
        pre = pre[:-2] + '│ '
    elif pre[-2:] == '└─':
        pre = pre[:-2] + '  '
    elif len(pre) < 2 or pre[-2] == '  ':
        pass
    else:
        raise Exception(f"Info structure failed: last two characters not recognized: '{pre}'")
    return pre


# printInfoStructure: print information from a list of Components in a structured way
#  @param label         Preceding text ending with '├─' or '└─'
#  @param compent_list  List of Components from which the info is printed
def print_components_info(pre, component_list):
    """
    This function accepts a list of Components and print their info in a structured tree-like format as follows:
    ├─Component0.print_components_info
    │ XXX
    └─Component1.print_components_info
      XXX
    """
    pre = update_pre(pre)
    for component in component_list[:-1]:
        component.print_components_info(pre + '├─')
    component_list[-1].print_components_info(pre + '└─')


# print box
def print_box(*args, layout=None, box_layout=None, **kwargs):
    """
    This functions prints adds a box around the string text
    :param args: str   The arguments to be printed
    :param layout:     The text layout to be used: plain, bold, underline, inverse, negative, warning, fail, grey, red,
                       green, yellow, flue, magenta, cyan, white, black
    :param box_layout: The box layout to be used: plain, bold, underline, inverse, negative, warning, fail, grey, red,
                       green, yellow, flue, magenta, cyan, white, black
    :param kwargs:   The kwargs of print_info
    :return: str
    """
    text = "".join(map(str, args))
    n = len(text)
    box_l = layout_style.get(box_layout)
    text_l = layout_style.get(layout)
    plain = layout_style.get('plain')
    text = text_l + text + box_l
    top = '\n' + box_l + '┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘' + plain
    print_info(top + mid + bottom, **kwargs)


# timer-function
@contextmanager
def quick_timer(name=None, t=0, n=0, ms=False):
    """
    Contextmanager that prints the time to execute a piece of code.

    Nothing is recorded/saved, this is meant mainly
    as a quick tool to check how long specific
    (parts of) calculations take during development.

    name: printed to identify the timer
    t, n: integers to specify number of tabs and new
          lines respectively
    ms:   print time in ms (default is s)


    example:
        with tools.quicktimer('test', ms=True):
            a = 1 / 3
    """
    start_time = time.time()
    yield
    elapsed_time = time.time() - start_time
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsed_time * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsed_time:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    print(s)


# initialization time measuring function
def time_initialize(initialize):
    def wrapper(*args):
        self = args[0]
        if not hasattr(self, 'init_time'):
            self.init_time = 0.0
        start_time = time.time()
        initialize(*args)
        self.init_time += time.time() - start_time

    return wrapper


# run time measuring function
def time_solve_solution_step(solve_solution_step):
    def wrapper(*args):
        self = args[0]
        if not hasattr(self, 'run_time'):
            self.run_time = 0.0
        start_time = time.time()
        interface_output = solve_solution_step(*args)
        self.run_time += time.time() - start_time
        return interface_output

    return wrapper


# save time measuring function
def time_save(output_solution_step):
    def wrapper(*args):
        self = args[0]
        if not hasattr(self, 'save_time'):
            self.save_time = 0.0
        start_time = time.time()
        interface_output = output_solution_step(*args)
        self.save_time += time.time() - start_time
        return interface_output

    return wrapper


# pass on parameters
def pass_on_parameters(settings_from, settings_to, keys):
    for key in keys:
        if key in settings_to:
            print_info(f'WARNING: parameter "{key}" is defined multiple times in JSON file', layout='warning')
        settings_to[key] = settings_from[key]


# compare bounding box of ModelParts
def check_bounding_box(mp_a, mp_b, directions, tol_center_warning=.02, tol_center_error=.1,
                       tol_minmax_warning=.1, tol_minmax_error=.3):
    """
    Use this function to compare the bounding boxes of 2 ModelParts.

    mp_a and mp_b are the two ModelParts
    and directions is a list of the mapping direction, for example ['x0', 'y0', 'z0'].

    There are 4 keyword arguments, to overwrite the
    default tolerances:
        tol_center_warning = 0.02
        tol_center_error = 0.1
        tol_minmax_warning = 0.1
        tol_minmax_error = 0.3

    Returns nothing.
    """

    for d in directions:
        if d not in ['x0', 'y0', 'z0']:
            raise ValueError(f'"{d}" is not a valid direction')
    if len(directions) > 3:
        raise ValueError(f'Too many directions given')
    if len(directions) == 0:
        raise ValueError(f'No directions specified')

    # extract coordinate data
    mp_a_coords = np.zeros(shape=(mp_a.size, len(directions)))
    mp_b_coords = np.zeros(shape=(mp_b.size, len(directions)))
    for i, d in enumerate(directions):
        mp_a_coords[:, i] = getattr(mp_a, d)
        mp_b_coords[:, i] = getattr(mp_b, d)

    # get bounding boxes
    mp_a_min = np.min(mp_a_coords, axis=0)
    mp_a_max = np.max(mp_a_coords, axis=0)
    mp_b_min = np.min(mp_b_coords, axis=0)
    mp_b_max = np.max(mp_b_coords, axis=0)
    mp_a_center = (mp_a_min + mp_a_max) / 2
    mp_b_center = (mp_b_min + mp_b_max) / 2

    # get reference distance (= average length of bounding box diagonal)
    mp_a_diag = mp_a_max - mp_a_min
    mp_b_diag = mp_b_max - mp_b_min
    d_ref = np.linalg.norm((mp_a_diag + mp_b_diag) / 2)

    # calculate errors on bounding boxes
    error_center = np.linalg.norm(mp_a_center - mp_b_center) / d_ref
    error_min = np.linalg.norm(mp_a_min - mp_b_min) / d_ref
    error_max = np.linalg.norm(mp_a_max - mp_b_max) / d_ref

    # raise warning or error if necessary
    msg_1 = f'ModelParts "{mp_a.name}", "{mp_b.name}": '
    msg_2 = ' values differ by '
    msg_3 = f'\n\t"{mp_a.name}": minimal values = {mp_a_min} and maximal values = {mp_a_max}' \
            f'\n\t"{mp_b.name}": minimal values = {mp_b_min} and maximal values = {mp_b_max}'

    msg = f'{msg_1}center{msg_2}{100 * error_center:.1f}%' + msg_3
    if error_center > tol_center_error:
        raise ValueError(msg)
    if error_center > tol_center_warning:
        warnings.warn(msg, Warning)

    msg = f'{msg_1}min{msg_2}{100 * error_min:.1f}%' + msg_3
    if error_min > tol_minmax_error:
        raise ValueError(msg)
    if error_min > tol_minmax_warning:
        warnings.warn(msg, Warning)

    msg = f'{msg_1}max{msg_2}{100 * error_max:.1f}%' + msg_3
    if error_max > tol_minmax_error:
        raise ValueError(msg)
    if error_max > tol_minmax_warning:
        warnings.warn(msg, Warning)


# import module from path
def import_module(module_name, path):
    """
    Loads module from a(n absolute) path
    :param module_name: string
    :param path: string
    :return: module
    """
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_env():
    """
    function to write all the environment variables as dict in a env.pickle file
    """
    with open('env.pickle', 'wb') as f:
        pickle.dump(dict(os.environ), f)


def get_solver_env(solver_module_name, working_dir):
    """
    Use this function to get all the environment variables corresponding to a solver_wrapper module.
    This uses a python-dict in the file coconut/solver_modules.py to get the module load command, which is then executed
    in a process. The process also runs coconut.tools.write_env() to store the resulting environment variables in a
    pickle file. Finally, pickle file is loaded and returned as a python-dict.

    @param solver_module_name: module name of the solver wrapper,
    e.g. coconut.coupling_components.solver_wrappers.fluent.v2023R1.
    e.g. fluent.v2023R1
    @type: str

    @param working_dir: working directory of the solver where the simulation is run.
    @type: str

    @return: environment variables as python-dict
    @rtype: dict
    """
    env_filename = 'env.pickle'
    pre_modules = 'coconut.coupling_components.solver_wrappers.'
    # remove pre_modules from the solver_module_name
    solver_name = solver_module_name.replace(pre_modules, '')

    # get the module load command for the solver
    solver_load_cmd = solver_modules.get_solver_cmd(solver_name)
    if not solver_load_cmd:
        solver_load_cmd = 'echo'

    # run the module load command and store the environment
    try:
        subprocess.check_call(
            f'{solver_load_cmd} && python3 -c "from coconut import tools;tools.write_env()"',
            shell=True, cwd=working_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        raise RuntimeError(f'Module load command for solver wrapper {solver_name} failed.')

    # load the environment variables and return as python-dict
    env_filepath = join(working_dir, env_filename)
    with open(env_filepath, 'rb') as f:
        env = pickle.load(f)
    os.remove(env_filepath)

    return env


def solver_available(solver_module_name):
    """
    @param solver_module_name: module name of the solver wrapper,
    e.g. coconut.coupling_components.solver_wrappers.fluent.v2023R1
    e.g. fluent.v2023R1
    @type: str

    @return: presence of solver env
    @rtype: bool
    """
    pre_modules = 'coconut.coupling_components.solver_wrappers.'
    # remove pre_modules from the solver_module_name
    solver_name = solver_module_name.replace(pre_modules, '')

    try:
        # get the module load command for the solver
        solver_modules.get_solver_cmd(solver_name)
    except KeyError:
        return False
    return True


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.saved_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.saved_path)


# try several times to remove file or directory
def rm_timed(path: str, sleep: float = 0.5, attempts: int = 100) -> None:
    """
    Tries to remove path several times
    :param str path: path to be removed
    :param float sleep: wait time after each attempt
    :param int attempts: number of attempts
    :return: None
    """
    for i in range(attempts):
        try:
            shutil.rmtree(path)
            break
        except OSError:
            time.sleep(sleep)
    if os.path.exists(path):
        print_info(f'Timed out removing {path}', layout='warning')
        shutil.rmtree(path)


# remove a key in a nested dictionary/list
def remove_recursively(key_to_remove, structure):
    if isinstance(structure, dict):
        structure.pop(key_to_remove, None)
        structure = structure.values()
    elif not isinstance(structure, list):
        return
    for value in structure:  # structure is either dict.values() or list
        remove_recursively(key_to_remove, value)

def flatten_concatenation(matrix):
    flat_list = []
    for row in matrix:
        flat_list += row
    return flat_list