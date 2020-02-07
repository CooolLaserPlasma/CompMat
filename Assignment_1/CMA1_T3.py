"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 3

Finished! :) kommentera klart bara
latest update: 2020-01-31
"""
import numpy as np
from numpy import pi, sqrt, exp
import matplotlib.pyplot as ppl
from scipy.linalg import eigh

#%%

def func(r):            # function containing diffrent potential contributions
    return -1/r   

def wavefunc_H(r, h):   # ground state wave function of hydrogen
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
A = M(N-1, h, func(ri))

w, v = eigh(A, eigvals=(0,0))
u = v.reshape(N-1)

n_s = u**2/(4*pi)
phi = sqrt(n_s)/ri
norm = sqrt(4*pi*phi.dot(phi)*h) 
phi = phi/norm

phi_analytic = wavefunc_H(ri, h)


ppl.figure()
ppl.plot(ri, phi, label='numerical')
ppl.plot(ri, phi_analytic, '--', label='analytic')
ppl.legend()
ppl.show()