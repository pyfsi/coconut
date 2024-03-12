from coconut.coupling_components.solver_wrappers.kratos_structure.kratos_structure import SolverWrapperKratosStructure
from coconut import tools

import json
import os
import numpy as np
import re


def create(parameters):
    return SolverWrapperKratosStructure91(parameters)


class SolverWrapperKratosStructure91(SolverWrapperKratosStructure):
    version = '91'

    def __init__(self, parameters):
        super().__init__(parameters)
        self.env = tools.get_solver_env(__name__, self.working_directory)
        self.check_software()
