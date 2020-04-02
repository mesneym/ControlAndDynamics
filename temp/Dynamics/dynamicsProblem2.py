
import numpy as np 

def transZ(n):
    a = np.identity(4)
    a[2,3]=n
    return a

def transX(n):
    a=np.identity(4)
    a[0,3]=n
    return a

def transY(n):
    a=np.identity(4)
    a[1,3]=n
    return a

def rotX(n):
    x=np.deg2rad(n)
    a=np.array([[1,0,0,0],
                [0,np.cos(x),-np.sin(x),0],
                [0,np.sin(x),np.cos(x),0 ],
                [0, 0,0,1]
               ])
    return a

def rotY(n):
    x=np.deg2rad(n)
    a=np.array([[np.cos(x),0,np.sin(x),0],
                [0,1,0,0],
                [-np.sin(x),0,np.cos(x),0 ],
                [0, 0,0,1]
               ])
    return a

def rotZ(n):
    x=np.deg2rad(n)
    a=np.array([[np.cos(x),-np.sin(x),0,0],
                [np.sin(x),np.cos(x),0,0],
                [0,0,1,0 ],
                [0,0,0,1]
               ])
    return a


def A_dh(theta,d,a,alpha):
   return  np.dot(rotZ(theta),
           np.dot(transZ(d),
           np.dot(transX(a),
           rotX(alpha))))




print('##############################################################')
print('                      Question 2a                              ')
print('##############################################################')
print(' ')

# q = np.array([90, -30, 60, 2])
q= np.array([150, 45 , 30 , 2])
l1 = 1.5
l2 = 1
l3 = 0.5
A1=A_dh(q[0]-90,l1,0,-30)
A2=A_dh(q[1]-90,l2,0,90)
A3=A_dh(q[2]+60,0,0,-90)
AN=A_dh(0,q[3]+l3,0,0)

T0_1 = A1
T0_2 = np.dot(T0_1,A2)
T0_3 = np.dot(T0_2,A3)
T0_n = np.dot(T0_3,AN)


z0 = np.array([0, 0, 1, 0])
z1 = np.dot(T0_1,z0)[:3]
z2 = np.dot(T0_2,z0)[:3]
z3 = np.dot(T0_3,z0)[:3]
z0 = z0[:3]

O0 = np.array([0, 0, 0, 1])
O1 = np.dot(T0_1,O0)[:3]
O2 = np.dot(T0_2,O0)[:3]
On = np.dot(T0_n,O0)[:3]
O0 = O0[:3]



j1 = np.concatenate((np.cross(z0,(On - O0)),z0))
j2 = np.concatenate((np.cross(z1,(On - O1)),z1))
j3 = np.concatenate((np.cross(z2,(On - O2)),z2))
j4 = np.concatenate((z3,O0))




J = np.column_stack((j1,j2,j3,j4))
F = np.array([1, 0, 0, 0, 0, 0])
T = np.round(np.dot(J.T,F),3)

print('Joint torques ')
print(T)


print('##############################################################')
print('                      Question 2b                              ')
print('##############################################################')
print(' ')

T= np.array([2, 3, -1, 0.5])
J= J[0:3,:]
pseudoInv = np.dot(np.linalg.inv(np.dot(J,J.T)),J)
print(np.dot(pseudoInv,T))



