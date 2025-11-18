from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

file = 'reactions_autocat.txt'
reactions = read_file(file_name= file)
k = [1]*8
initial = [1000]*8
V = 1000

n_iterations = 1000
method = "Protocell"

abundances, times, V = chemistry(method= method, iterations= n_iterations,
                                 reactions= reactions, initial_food= initial,
                                 k= k, V= V)

data_file = open('example','ab')
to_save = np.vstack((abundances.T,times,V))
pickle.dump(to_save, data_file)
data_file.close()

print(to_save.shape)
print(to_save.T.shape)
print(to_save.T[:100,-2])