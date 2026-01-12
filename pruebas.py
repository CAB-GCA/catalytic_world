import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('..')  # allow importing from parent directory
from fun_gilles import *
import pickle
import pandas as pd

a, t, v = chemistry(method= "Protocell",
                    iterations=1000,
                    file= "examples/reactions_autocat.txt",
                    initial_food=[1.e+02, 1.e+02, 1.e+02, 1.e+02, 1.e-02, 0.e+00, 0.e+00, 0.e+00],
                    k= [1] + [1e-4] + [1]*6,
                    V= 100)
print(np.sum(a[:,4]==0))
