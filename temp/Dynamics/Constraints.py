import numpy as np

print('##############################################################')
print('                      Question 1a                              ')
print('##############################################################')
print(' ')

a = [[1, 0, 0, 0],
     [0, 2, 3, 0],
     [1, 0, 0, 5]]
J = np.array(a)
print(np.linalg.inv(np.dot(J,J.T)))


print('##############################################################')
print('                      Question 1c                              ')
print('##############################################################')
print(' ')

T1 = np.array([2, 4, 6, 5])
T2 = np.array([2, 5, 6, 5])
pseudoInv = np.dot(np.linalg.inv(np.dot(J,J.T)),J)
print(np.dot(pseudoInv,T1))
print(np.dot(pseudoInv,T2)) 


