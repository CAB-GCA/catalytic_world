import numpy as np
from fun_gilles import *


reactions = read_file('reactions_prueba.txt')
m_reactions = c_matrix(reactions, obtain_species(reactions))
print(m_reactions)
print(update_a(m_reactions))

reactions2 =  read_file('reactions_0.txt')
m_reactions2 = c_matrix(reactions2, obtain_species(reactions2))
print(m_reactions2)
print(update_a(m_reactions2))

