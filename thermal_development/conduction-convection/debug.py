import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle

with open('case_results.pickle', 'rb') as file:
    results = pickle.load(file)

print(results)