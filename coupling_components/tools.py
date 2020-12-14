import time
from contextlib import contextmanager


def create_instance(settings):
    # create instance of given class based on settings dict
    object_type = settings['type']
    object_module = __import__('coconut.coupling_components.' + object_type, fromlist=[object_type])
    return object_module.create(settings)


# PrintInfo: Printing information with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintInfo(label, *args):
    print(label, " ".join(map(str, args)))


# PrintWarning: Printing a warning with a label
#
#  @param label         The label for the print
#  @param args          The arguments to be printed
def PrintWarning(label, *args):
    print(label, " ".join(map(str, args)))


# Class contains definition of colors. This is to be used as a struct.
#
# Example usage: print(bcolors.HEADER + "This is a header in header color" + bcolors.ENDC)
# IMPORTANT: The end of the print statement should always contain bcolor.ENDC
class bcolors:
    HEADER    = '\033[95m'
    BLUE      = '\033[94m'
    GREEN     = '\033[92m'
    MEGENTA   = '\033[96m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


class LayoutStyles:
    styles = {'header':    '\033[95m',
              'blue':      '\033[94m',
              'green':     '\033[92m',
              'red':       '\033[96m',
              'warning':   '\033[1;33;41m',
              'fail':      '\033[91m',
              'bold':      '\033[1m',
              'underline': '\033[4m',
              'plain':     '\033[0m'
              }

    def get(self, layout):
        self.check(layout)
        return self.styles[layout]

    def check(self, layout):
        if layout not in self.styles:
            raise ValueError("Layout style is not implemented, correct layout styles are:"
                             "header, blue, green, red, warning, fail, bold, underline and plain.")


layout_style = LayoutStyles()


# Print: Printing with color
#
#  @param args          The arguments to be printed
#  @param layout        The layout to be used: header, blue, green, red, warning, fail, bold, underline or plain
def Print(*args, layout=None):
    if layout is None:
        print("".join(map(str, args)))
    else:
        print(layout_style.get(layout), "".join(map(str, args)), layout_style.get('plain'))


# UpdatePre: Update preceding text, used in structure printing
#  @param pre         Preceding text ending with '├─' or '└─'
def UpdatePre(pre):
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


# PrintInfoStructure: Print information from a list of Comonents in a strucutred way
#  @param label         Preceding text ending with '├─' or '└─'
#  @param compent_list  List of Components from which the info is printed
def PrintComponentsInfo(pre, component_list):
    """
    This function accepts a list of Components and print their info in a structured tree-like format as such:
    ├─Component0.PrintInfo
    │ XXX
    └─Component1.PrintInfo
      XXX
    """
    pre = UpdatePre(pre)
    for component in component_list[:-1]:
        component.PrintComponentsInfo(pre + '├─')
    component_list[-1].PrintComponentsInfo(pre + '└─')


# Timer-function
@contextmanager
def quicktimer(name=None, t=0, n=0, ms=False):
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
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    if ms:
        s = '\n' * n + '\t' * t + f'{elapsedTime * 1000:.2f}ms'
        s.replace(',', ' ')
    else:
        s = '\n' * n + '\t' * t + f'{elapsedTime:.1f}s'
    if name is not None:
        s += f' - {name}'
    s += '\n' * n
    Print(s)


# Run time measuring function
def TimeSolveSolutionStep(SolveSolutionStep):
    def wrapper(*args):
        self = args[0]
        if not hasattr(self, 'run_time'):
            self.run_time = 0.0
        start_time = time.time()
        interface_output = SolveSolutionStep(*args)
        self.run_time += time.time() - start_time
        return interface_output

    return wrapper
