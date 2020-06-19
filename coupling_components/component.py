from coconut.coupling_components import tools


class Component(object):
    def __init__(self):
        self.initialized = False
        self.initialized_solution_step = False

    def Initialize(self):
        if self.initialized:
            raise Exception("Already initialized")
        else:
            self.initialized = True

    def Finalize(self):
        if self.initialized_solution_step:
            raise Exception("Solution step not finalized")
        if self.initialized:
            self.initialized = False
        else:
            raise Exception("Not initialized")

    def InitializeSolutionStep(self):
        if self.initialized:
            if self.initialized_solution_step:
                raise Exception("Already solution step initialized")
            else:
                self.initialized_solution_step = True
        else:
            raise Exception("Not initialized")

    def FinalizeSolutionStep(self):
        if self.initialized:
            if self.initialized_solution_step:
                self.initialized_solution_step = False
            else:
                raise Exception("Solution step not initialized")
        else:
            raise Exception("Not initialized")

    def OutputSolutionStep(self):
        pass

    def Check(self):
        pass

    def PrintInfo(self, pre):
        tools.Print(pre, "The component ", self.__class__.__name__)