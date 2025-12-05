from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

file = 'examples/reactions_autocat.txt'
reactions = read_file(file_name= file)

n_iterations = 2e5
method = "Protocell" # Gillespie or Deterministic
# Reaction constants:
k = [1]*8 # len(k)= # de reacciones
# Volume:
V = 1000
initial_food = [1000]*4 + [0]*4 # initial molecules number
final_volume = []
k_change = np.logspace(-3,4,7)

for k_i in k_change:
    k[0] = k_i
    print(f"Performing simulation for k = {k_i}")
    abundances, times, volumes = chemistry(method, n_iterations, file, initial_food, k, V)
    print(f"Simulation ended. Final volume = {volumes[-1]}")
    
    try:
        print(f"Volume at t= 10 = {volumes[times>10][0]}")
        final_volume.append(volumes[times>10][0])
    except IndexError:
        final_volume.append(volumes[-1])

# # data_file = open('example','ab')
# # to_save = np.vstack((abundances.T,times,V))
# # pickle.dump(to_save, data_file)
# # data_file.close()
# # print(to_save.shape)
# # print(to_save.T.shape)
# # print(to_save.T[:100,-2])

# def block_statistics(concentration, W):
#     """
#     Calculates the mean and standard deviation for non-overlapping blocks.
#     W is the block size (e.g., 100).
#     """
#     N = len(concentration)
    
#     # 1. Determine the number of complete blocks
#     num_complete_blocks = N // W
    
#     # 2. Slice the data to discard the remainder (if any)
#     sliced_concentration = concentration[:num_complete_blocks * W]

#     # 3. Reshape the data into blocks (e.g., 100 blocks of 100 iterations)
#     # Each row is now a block of W iterations.
#     blocked_data = sliced_concentration.reshape(num_complete_blocks, W)

#     # 4. Calculate the mean and standard deviation along axis=1 (across the W iterations)
#     block_mean = np.mean(blocked_data, axis=1)
#     block_std = np.std(blocked_data, axis=1)
    
#     return block_mean, block_std

# mean, std = (block_statistics(abundances[-1000:,0], 50))
# print(mean)
# print(std)


# from numpy.linalg import matrix_rank


# matrix_rank(np.eye(4)) # Full rank matrix
# print(np.eye(4))
# print(matrix_rank(np.eye(4)))
# I=np.eye(4); I[-1,-1] = 0. # rank deficient matrix
# print(matrix_rank(I))
# print(matrix_rank(np.ones((4,))) )# 1 dimension - rank 1 unless all 0
# print(matrix_rank(np.zeros((4,))))

# print(reactions)
# print(reactions[reactions[:,-1]!='4'])
# c_matrix_wo_food = c_matrix(reactions[reactions[:,-1]!='4'], obtain_species(reactions))
# print(c_matrix_wo_food)
# print(matrix_rank(c_matrix_wo_food))

# print(check_thermodynamics(reactions))


# reactions = read_file('examples/incompatible.txt')
# print(check_thermodynamics(reactions))
# print(1000%102)


# def block_statistics(concentration):
#     """
#     Calculates the mean and standard deviation for non-overlapping blocks.
#     W is the block size (e.g., 100).
#     """
#     block_mean = np.mean(concentration, axis=0)
#     block_std = np.std(concentration, axis=0)
    
#     return block_mean, block_std

# print(block_statistics(abundances[:10]))

# print(sum([True, False, True]))