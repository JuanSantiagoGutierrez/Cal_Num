import numpy as np


A = np.array([
    [1,0,0,0,0],
    [1,1000,0,0,0],
    [1,2000,2e6,0,0],
    [1,3000,6e6,6e9,0],
    [1,4000,12e6,2.4e10,2.4e13]
    ])
b = np.array([175,350,330,300,150])

coeff_X = np.linalg.solve(A, b)

print("Coeff: ", coeff_X)

list_x = [1000, 2000, 3000, 4000, 5000]
n = 1
suma = 0
for i in range(len(coeff_X)):
    suma += n * coeff_X[i]
    n = n * (1500 - list_x[i])

print(suma)
