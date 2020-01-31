"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 4
"""
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as ppl
from scipy.linalg import eigh, inv

#%%
# From Task 2
def task2(n_s, r, h, N):
    f = -4*pi*r*n_s  
    M = -2*np.identity(N-1) + np.diag(np.ones(N-2),1)+np.diag(np.ones(N-2),-1)
    U_0 = h**2*inv(M).dot(f)  
    return (U_0 + (r/r[-1]))/r

# from Task 3
def task3(f, r, h, N):
    
    def M(n,h,v):
        out = -2*np.identity(n)+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
        out = -out/(2*h**2) + np.diag(v)
        return out
    
    A = M(N-1, h, f)

    w, v = eigh(A, eigvals=(0,0))
    u = v.reshape(N-1)

    n_s = u**2/(4*pi)
    phi = sqrt(n_s)/r
    norm = sqrt(4*pi*phi.dot(phi)*h) 
    phi = phi/norm
    return (w, phi)

# ground state electron density for hydrogen atom
def n_g_H(r, h):
    n = np.exp(-2*r)                   
    norm = 4*pi*h*(r**2).dot(n)
    return n/norm

#%%
N = 10
r = np.linspace(0, 10, N)
ri = r[1:]
h = r[1]-r[0]

# guess some density
n_s = n_g_H(ri, h) 
# Calculate initial values
phi = sqrt(n_s)
f = -1/ri + task2(n_s, ri, h, N)
E_0, phi = task3(f, ri, h, N)  

#%%
eV = 1/27.211     # Hartree
E_dif_tol = 10**(-5)*eV  #10**(-5)*eV

i_max = 20

for i in range(i_max):
    # solve Poisson (as in task2) n_s --> U --> V_sH
    V_sH = task2(n_s, ri, h, N)
    # create input function for task3 as
    f = -1/ri + V_sH
    # solve eigenv.p as in task 3 V_sH=V_H --> u, e
    E_0_temp, phi = task3(f, ri, h, N)                    
    n_s = abs(phi)**2
    # compare energy with previous energy  
    E_dif = abs(E_0 - E_0_temp)
    E_0 = E_0_temp
    if E_dif < E_dif_tol:
        print('E_dif_tol reached')
        print('Number of iterations: ' + str(i))
        print('Ground state energy: ' + str(E_0))
        break


#%%
        

ppl.figure()
ppl.plot(ri, phi, label='numerical')
ppl.plot(ri, phi2, '--', label='from task3 func')
ppl.legend()
ppl.show()