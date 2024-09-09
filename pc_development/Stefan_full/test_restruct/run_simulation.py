import coconut

import json


# import parameters
parameter_file_name = "parameters_test.json"
with open(parameter_file_name, 'r') as parameter_file:
    parameters = json.load(parameter_file)

simulation = coconut.Analysis(parameters)
simulation.run()
