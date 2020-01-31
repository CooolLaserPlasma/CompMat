"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 1
testfile
"""
import numpy as np
from scipy.linalg import eigh, inv
from numpy import pi, sqrt, exp
import matplotlib.pyplot as ppl

#%% 1
a_1 = 0.298073  # 0.297104   #α 1 = 0.298 073
a_2 = 1.242567  # 1.236745   #α 2 = 1.242 567
a_3 = 5.782948  # 5.749982   #α 3 = 5.782 948 
a_4 = 38.474970 # 38.216677  #α 4 = 38.474 970.

# function for genarating elements of the Q-matrix (OK! 2020-01-29)
def Q_prqs(a_p, a_r, a_q, a_s):
    return 2*pi**(5/2)/((a_p+a_q)*(a_r+a_s)*sqrt(a_p+a_q+a_r+a_s))



#%% 2
# generating the Q-matrix (OK! 2020-01-29)
a = np.array([a_1, a_2, a_3, a_4])
Q = np.zeros([4, 4, 4, 4])



for p in range(4):
    for r in range(4):
        for q in range(4):
            for s in range(4):
                Q[p, r, q, s] = Q_prqs(a[p], a[r], a[q], a[s])

#%% 3
# function for generating elements of the S-matrix (OK! 2020-01-29)
def S_pq(a_p, a_q):
    return (pi/(a_p + a_q))**(3/2)



# Generate the S-matrix (OK! 2020-01-29)
S = np.zeros([4, 4])

for p in range(4):
    for q in range(4):
        S[p,q] = S_pq(a[p], a[q])


#%% 4
# function for generating the elements of the h-matrix
def h_pq(a_p, a_q):
    alpha = a_p + a_q
    t_1 = 3*a_q*a_p*pi**(3/2)/(alpha**(5/3))
    t_2 = 4*pi/alpha
    return t_1 - t_2

# Generate h-matrix
h = np.zeros([4, 4])

for p in range(4):
    for q in range(4):
        h[p,q] = h_pq(a[p], a[q])
        
#%% 5
# initial values for the coefficients C_p
C_1 = 1
C_2 = 1
C_3 = 7
C_4 = 1

C = np.array([C_1, C_2, C_3, C_4])

# Normalize C 
norm = C.dot(S.dot(C))
C = C/norm

# create F-matrix 
F = h + C.dot(Q.dot(C))

#%%
# solve generelized eigenvalue equation
w, v = eigh(F, S, eigvals=(0,0))
C = v.reshape(4)

#%%
# 
eV = 1/27.211     # Hartree
E_dif_tol = 10**(-5)*eV  #10**(-5)*eV
# Calculate ground state energy
E_g = 2*C.dot(h.dot(C)) + C.dot(C.dot(C.dot(C.dot(Q))))

#%%
i_max = 20

for i in range(i_max):
    # normalize C
    norm = C.dot(S.dot(C))
    C = C/norm
    # beräkna F
    F = h + C.dot(Q.dot(C))
    # solve eigh
    w, v = eigh(F, S, eigvals=(0,0))
    C = v.reshape(4)
    # Calculate new energy  
    E_g_temp = 2*C.dot(h.dot(C)) + C.dot(C.dot(C.dot(C.dot(Q))))
    # compare energy with previous energy  
    E_dif = abs(E_g - E_g_temp)
    #print(E_dif)
    E_g = E_g_temp
    print(E_g)
    # if dif < dif_tol
    #       exit
    if E_dif < E_dif_tol:
        print('E_dif_tol reached')
        print('Number of iterations: ' + str(i))
        print('Ground state energy: ' + str(E_g))
        break


##
#E_goal = −2.855 160 38
#%%
# create wave-function

def phi(x, C, a):
    return C[0]*exp(-a[0]*x**2) + C[1]*exp(-a[1]*x**2) + C[2]*exp(-a[2]*x**2) + C[3]*exp(-a[3]*x**2)

X = np.linspace(0, 5, 100)
Y = phi(X, C, a)

ppl.plot(X, Y) 

#def phi(x, C, a):
#    return C[0]*exp**(-a[0]*x.**2)
    
        
    
#%% 
"""
###############################################################################
Test
###############################################################################
"""
# 2020-01-28
A = np.zeros([2, 2, 2, 2])
p = 0

for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                p = p + 1
                A[i, j, k, l] = p

#%% This works!

A_temp = np.zeros([2, 2])

for r in range(2):
    for s in range(2):
        A_temp = A_temp + A[:, r, :, s]



B = A[:, 0, :, 0] + A[:, 0, :, 1] + A[:, 1, :, 0] + A[:, 1, :, 1]

C = sum(A[:, p, :, s] for p in range(2) for s in range(2))   # This works aswell!!:) 

#%%
C = np.array([1, 2])
B = A[:, 0, :, 0]*C[0]*C[0] + A[:, 0, :, 1]*C[0]*C[1] + A[:, 1, :, 0]*C[1]*C[0] + A[:, 1, :, 1]*C[1]*C[1]

D = C.dot(A.dot(C))

#%%
# Solving eigenvalue problems

V = np.array([[1, 5], [3, 2]])
S = inv(V.T).dot(inv(V))

H_p = np.array([[2, 1], [1, 2]])
C_1p = np.array([1, -1])
C_2p = np.array([1, 1])
E_1 = 1
E_2 = 3
C_1 = V.dot(C_1p)
C_2 = V.dot(C_2p)
H = inv(V.T).dot(H_p.dot(inv(V)))


w, v = eigh(H, S)
w_2, v_2 = eigh(H, S, eigvals=(0,0))

# OBS skalfaktor!
# C_1 = v[:,0]
# C_2 = v[:,1]

#%%
# test looop
E_g = 10
i_max = 20
E_dif_tol = 0.2

for i in range(i_max):
    d = 1/(1+i)
    #print(d)
    E_g_temp = E_g+d       #calculate new energy
    E_dif = abs(E_g - E_g_temp)
    #print(E_dif)
    if E_dif < E_dif_tol:
        print('E_dif_tol reached')
        break
    else:
        E_g = E_g_temp
        #print(E_g)

# while loop
dif_tol = 2
E_0 = 20
E = 0
dif = E_0-E

#%%
while dif > dif_tol:
    dif = abs(E_0-E)
    E = E + 1
    print(dif)
    #E = E_temp
    #print(E)

#%%
"""
LOG
-----------------------------------------------
2020-01-28:
cells (1-3) all OK!

2020-01-29
did not save yesterday but recovered stuff from history.py file

F and S are symetric OK!

"""
#%%
"""
(12*pi*a_1*x^2-8*pi*a_1^2*x^4-2/x)e^(-(a_1+a_2)x^2)
"""