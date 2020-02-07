"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 1

Finished! :)
latest update: 2020-01-31
"""
import numpy as np
from scipy.linalg import eigh
from numpy import pi, sqrt, exp
import matplotlib.pyplot as ppl

#%%
# function for genarating elements of the Q-matrix 
def Q_prqs(a_p, a_r, a_q, a_s):
    return 2*pi**(5/2)/((a_p+a_q)*(a_r+a_s)*sqrt(a_p+a_q+a_r+a_s))

# function for generating elements of the S-matrix 
def S_pq(a_p, a_q):
    return (pi/(a_p + a_q))**(3/2)

# function for generating the elements of the h-matrix
def h_pq(a_p, a_q):
    alpha = a_p + a_q
    t_1 = 3*a_q*a_p*pi**(3/2)/(alpha**(5/2))
    t_2 = 4*pi/alpha
    return t_1 - t_2

#%%
a_1 = 0.297104   
a_2 = 1.236745   
a_3 = 5.749982   
a_4 = 38.216677  

a = np.array([a_1, a_2, a_3, a_4])

# generating the Q-matrix 
Q = np.zeros([4, 4, 4, 4])

for p in range(4):
    for r in range(4):
        for q in range(4):
            for s in range(4):
                Q[p, r, q, s] = Q_prqs(a[p], a[r], a[q], a[s])

# Generate the S-matrix 
S = np.zeros([4, 4])

for p in range(4):
    for q in range(4):
        S[p,q] = S_pq(a[p], a[q])

# Generate h-matrix
h = np.zeros([4, 4])

for p in range(4):
    for q in range(4):
        h[p,q] = h_pq(a[p], a[q])

#%%       
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

# solve generelized eigenvalue equation
w, v = eigh(F, S, eigvals=(0,0))
C = v.reshape(4)

eV = 1/27.211                                                       # Hartree
E_dif_tol = 10**(-5)*eV  
E_g = 2*C.dot(h.dot(C)) + C.dot(C.dot(C.dot(C.dot(Q))))             # Calculate ground state energy

i_max = 20

for i in range(i_max):
    norm = C.dot(S.dot(C))                                          # normalize C
    C = C/norm
    F = h + C.dot(Q.dot(C))                                         # soolve eigenvalue problem
    w, v = eigh(F, S, eigvals=(0,0))
    C = v.reshape(4)                                                #  Calculate new energy  
    E_g_temp = 2*C.dot(h.dot(C)) + C.dot(C.dot(C.dot(C.dot(Q)))) 
    E_dif = abs(E_g - E_g_temp)                                     # compare energy with previous energy 
    E_g = E_g_temp
    if E_dif < E_dif_tol:
        print('E_dif_tol reached')
        print('Number of iterations: ' + str(i))
        print('Ground state energy: ' + str(E_g))
        break

#%%
# create wave-function
def phi(x, C, a):
    return C[0]*exp(-a[0]*x**2) + C[1]*exp(-a[1]*x**2) + C[2]*exp(-a[2]*x**2) + C[3]*exp(-a[3]*x**2)

X = np.linspace(0, 5, 100)
Y = phi(X, C, a)

fig, ax = ppl.subplots()
ax.plot(X, Y)
ax.set_ylabel('Wavefunction')
ax.set_xlabel('r')
    
        