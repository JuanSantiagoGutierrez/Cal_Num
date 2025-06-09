# Tarea ECN - Problema aplicado autovalores

import numpy as np

E = 10e9 # [Pa]
i = 1.25e-5 # [m⁴]
L = 3 # [m]

A = np.array((
    [-2, 1, 0, 0, 0],
    [1, -2, 1, 0, 0],
    [0, 1, -2, 1, 0],
    [0, 0, 1, -2, 1],
    [0, 0, 0, 1, -2]))

eigA , veceigA = np.linalg.eig(A) 

k = - np.pow(L, 2) / (np.pow(E * i, 2) * 36)

P_list = []
for i in range(5):
    P = np.sqrt((eigA[i]/k))
    P_list.append(P)
    
print(min(P_list))