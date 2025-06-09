import numpy as np
A = np.array([[3,0,0,0],[1,2,0,0],[3,4,1,0],[1,2,3,4]])

b = np.array([1,0,0,0])

def sustProg(A:np.matrix, b:np.array):
    X = np.zeros((len(b)))
    X[0] = b[0]/A[0][0]
    X[1] = b[1] - A[1][0]*X[1]
    for i in range(2,len(b)):
        X[i] = (b[i] -  X[0:i].dot(A[i][0:i]))
    return X



B = np.matrix(A)
#print(A.diagonal())
#print(B.diagonal())
