import numpy as np

P = [[1,0,0],
     [0,0,1],
     [0,1,0]]


A = [[2,3,4],
     [5,6,7],
     [8,9,10]]

P = np.array(P)
A = np.array(A)

A_local = P.dot(A).dot(P)
A_global = P.dot(A_local).dot(P)

if A_global.any() == A.any():
    print('Ok')
else:
    print('Not Ok')

