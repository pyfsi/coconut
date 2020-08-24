from coconut.data_structure import KratosUnittest
from coconut.tests.solver_wrappers.fluent.test_v2019R1 import TestSolverWrapperFluent2019R1


def print_box(text):
    n = len(text)
    top = '\n┌─' + n * '─' + '─┐'
    mid = '\n│ ' + text + ' │'
    bottom = '\n└─' + n * '─' + '─┘'
    print(top + mid + bottom)

# *** TO DO: put two test cases in 1 folder to reduce number of folders

class TestSolverWrapperFluent2019R2(TestSolverWrapperFluent2019R1):
    # warning: because of inheritance, this class cannot be tested independently

    version = '2019R2'

    def test_solver_wrapper_fluent_2019R2(self):
        super().test_solver_wrapper_fluent_2019R1_tube2d()
        super().test_solver_wrapper_fluent_2019R1_tube3d()


if __name__ == '__main__':
    KratosUnittest.main()
