"""
TIF320 - Computational materials and molecular physics
Assignment 1
Task 2

Finished! :)
latest update: 2020-01-31
"""
import numpy as np
from numpy import pi
import matplotlib.pyplot as ppl
from scipy.linalg import inv

#%%

def n_s(r, h):
    n = np.exp(-2*r)                   # ground state electron density for hydrogen atom
    norm = 4*pi*h*(r**2).dot(n)
    return n/norm

a = 0 
b = 10
N = 1000
r = np.linspace(a,b,N)
ri = r[1:]
h = r[1]-r[0]
f = -4*pi*ri*n_s(ri, h)  


M = -2*np.identity(N-1) + np.diag(np.ones(N-2),1)+np.diag(np.ones(N-2),-1)

U_0 = h**2*inv(M).dot(f)
U = U_0 + (ri/ri[-1])   

V_analytic = 1/ri - (1+1/ri)*np.exp(-2*ri)    # Hartree potential

ppl.figure()
ppl.plot(ri, U/ri, label='numerical')
ppl.plot(ri, V_analytic, '--', label='analytic')
ppl.legend()
ppl.show()