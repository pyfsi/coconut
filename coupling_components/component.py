from coconut import tools


class Component(object):
    mapped = False
    dummy = False

    def __init__(self):
        self.initialized = False
        self.initialized_solution_step = False

    def initialize(self, *args):
        if self.initialized:
            raise Exception("Already initialized")
        else:
            self.initialized = True

    def finalize(self):
        if self.initialized_solution_step:
            raise Exception("Solution step not finalized")
        if self.initialized:
            self.initialized = False
        else:
            raise Exception("Not initialized")

    def initialize_solution_step(self):
        if self.initialized:
            if self.initialized_solution_step:
                raise Exception("Already solution step initialized")
            else:
                self.initialized_solution_step = True
        else:
            raise Exception("Not initialized")

    def finalize_solution_step(self):
        if self.initialized:
            if self.initialized_solution_step:
                self.initialized_solution_step = False
            else:
                raise Exception("Solution step not initialized")
        else:
            raise Exception("Not initialized")

    def output_solution_step(self):
        pass

    def print_components_info(self, pre):
        tools.print_info(pre, "The component ", self.__class__.__name__)
