from fun_gilles import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

file = 'examples/reactions_autocat.txt'
reactions = read_file(file_name= file)
k = [1]*8
initial = [1000]*8
V = 1000

# n_iterations = 1000
# method = "Protocell"

# abundances, times, V = chemistry(method= method, iterations= n_iterations,
#                                  reactions= reactions, initial_food= initial,
#                                  k= k, V= V)

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


from numpy.linalg import matrix_rank


matrix_rank(np.eye(4)) # Full rank matrix
print(np.eye(4))
print(matrix_rank(np.eye(4)))
I=np.eye(4); I[-1,-1] = 0. # rank deficient matrix
print(matrix_rank(I))
print(matrix_rank(np.ones((4,))) )# 1 dimension - rank 1 unless all 0
print(matrix_rank(np.zeros((4,))))

print(reactions)
print(reactions[reactions[:,-1]!='4'])
c_matrix_wo_food = c_matrix(reactions[reactions[:,-1]!='4'], obtain_species(reactions))
print(c_matrix_wo_food)
print(matrix_rank(c_matrix_wo_food))

print(check_thermodynamics(reactions))


reactions = read_file('examples/incompatible.txt')
print(check_thermodynamics(reactions))
abundances, V, times = chemistry("Protocell", 100, reactions, [10,10,10], [1]*6, 100)