"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 3

testfile 
latest update: 2020-01-31
"""
import numpy as np
from numpy import pi, sqrt, exp
import matplotlib.pyplot as ppl
from scipy.linalg import eigh

#%%

def func(r):
    return -1/r   

def wavefunc_H(r, h):
    wavef = exp(-r)
    norm = sqrt(4*pi*wavef.dot(wavef)*h) 
    return wavef/norm

def M(n,h,v):
    out = -2*np.identity(n)+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
    out = -out/(2*h**2) + np.diag(v)
    return out

N = 1000
r = np.linspace(0, 10, N)
ri = r[1:]
h = r[1]-r[0]
#M = (-2*np.identity(N-1) + np.diag(np.ones(N-2),1)+np.diag(np.ones(N-2),-1))/(2*h**2) 
#F = np.diag(f(ri))
A = M(N-1, h, func(ri)) #M+F

w, v = eigh(A, eigvals=(0,0))
u = v.reshape(N-1)

n_s = u**2/(4*pi)
phi = sqrt(n_s)/ri
norm = sqrt(4*pi*phi.dot(phi)*h) 
phi = phi/norm

phi_analytic = wavefunc_H(ri, h)

#E_0 = 2*w - 

ppl.figure()
ppl.plot(ri, phi, label='numerical')
ppl.plot(ri, phi_analytic, '--', label='analytic')
ppl.legend()
ppl.show()

#%%
"""
(s. 36 Thijssen)
Check 1 Fortunately, we again have an exact answer for the ground state energy:
this should be equal to −0.5 hartree = 13.6058 eV, and, if your program contains
no errors, you should find −0.499 278 hartree, which is amazingly good if you
realise that only four functions have been taken into account.
"""
E_gH = -0.5
#%% Making a function of this task
def func(r):
    return -1/r   

N = 1000
r = np.linspace(0, 10, N)
ri = r[1:]
h = r[1]-r[0]

f = func(ri)

"""
Functin OK!! :) 2020-01-31
"""
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

E_0, phi2 = task3(f, ri, h, N)

#phi2 = sqrt(n_s)/ri
#norm = sqrt(4*pi*phi2.dot(phi2)*h) 
#phi2 = phi2/norm

ppl.figure()
#ppl.plot(ri, phi, label='numerical')
ppl.plot(ri, phi2, '--', label='from task3 func')
ppl.legend()
ppl.show()