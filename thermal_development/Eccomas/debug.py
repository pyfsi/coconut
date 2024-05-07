import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

path_1 = 'case_results.pickle'
path_2 = 'test_single_solver/case_results.pickle'

with open(path_1, 'rb') as file:
    results = pickle.load(file)

print(results)