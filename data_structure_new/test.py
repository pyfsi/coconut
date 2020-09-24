from coconut.coupling_components.tools import quicktimer

from coconut.data_structure_new.model import Model
from coconut.data_structure_new.interface import Interface

from coconut.data_structure.Model import Model as ModelOld
from coconut.coupling_components.interface import Interface as InterfaceOld
from coconut.data_structure.Parameters import Parameters
from coconut import data_structure

import numpy as np
import json
import timeit


# >> testing code functionality
if 0:
    model = Model()
    for i in range(3):
        coords = np.random.rand(5 + i, 3)
        model.create_model_part(f'mp{i}', coords[:, 0], coords[:, 1], coords[:, 2])

    print(model)
    mp1 = model.get_model_part('mp1')

    print(mp1)

    par = {'mp0': ['pressure'],
           'mp1': ['pressure', 'traction']}  # *** order correct?? orderd dicts??

    interface = Interface(par, model)
    print(interface)

    print(interface.get_interface_data())

    interface.set_variable_data('mp1', 'pressure', np.random.rand(6, 1))

    print(interface.get_variable_data('mp1', 'pressure'))
    print(interface.get_interface_data())
    print('\n' * 5)


# >> compare speed of old and new data structure
if 0:
    n = 10000
    coords = np.random.rand(n, 3)

    par = {'mp0': ['pressure'],
           'mp1': ['pressure', 'traction']}
    tmp = {'mp0': ['PRESSURE'],
           'mp1': ['PRESSURE', 'TRACTION']}
    par_old = Parameters(json.dumps(tmp))

    # setup with new data structure
    model = Model()
    for model_part_name in par:
        model.create_model_part(model_part_name, *np.hsplit(coords, 3))
    interface = Interface(par, model)

    # setup with old data structure
    model_old = ModelOld()
    for model_part_name in par_old.keys():
        mp = model_old.CreateModelPart(model_part_name)
        for var in par_old[model_part_name].list():
            mp.AddNodalSolutionStepVariable(vars(data_structure)[var.GetString()])
        for i in range(n):
            mp.CreateNewNode(i, coords[i, 0], coords[i, 1], coords[i, 2])
    interface_old = InterfaceOld(model_old, par_old)

    print('\nget and set numpy array in interface')
    with quicktimer('new', t=1, ms=True):
        data = interface.get_interface_data()
        interface.set_interface_data(np.random.rand(*data.shape))
    with quicktimer('old', t=1, ms=True):
        data_old = interface_old.GetNumpyArray()
        interface_old.SetNumpyArray(np.random.rand(*data_old.shape))

    print('\ncopy interface')
    with quicktimer('new', t=1, ms=True):
        interface_2 = interface.copy()
    with quicktimer('old', t=1, ms=True):
        interface_old_2 = interface_old.deepcopy()

    print('\nadd interfaces')
    with quicktimer('new', t=1, ms=True):
        interface + interface_2
    with quicktimer('old', t=1, ms=True):
        interface_old + interface_old_2

# test speed of methods
if 1:
    setup = \
    """
from coconut.data_structure_new.model import Model
from coconut.data_structure_new.interface import Interface
import numpy as np

n = 1000
coords = np.random.rand(n, 3)

par = {'mp0': ['pressure'],
   'mp1': ['pressure', 'traction']}

model = Model()
for model_part_name in par:
    model.create_model_part(model_part_name, *np.hsplit(coords, 3))
interface = Interface(par, model)

data = interface.get_interface_data()
data = np.random.rand(*data.shape)

interface2 = interface.copy()
    """

    code1 = """interface.get_interface_data()"""
    code2 = """interface.set_interface_data(data)"""
    code3 = """interface.copy_old()"""
    code4 = """interface.copy()"""
    code5 = """interface + interface2"""
    code6 = """interface += interface2"""


    m = 10000
    codes = [code1, code2, code3, code4, code5, code6]
    # codes = [code3, code4]
    for code in codes:
        print(f'{code}:\n\t{timeit.timeit(setup=setup, stmt=code, number=m) / m * 1e6:} micro-s')

