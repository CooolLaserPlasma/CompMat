#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trying out ASE
"""
from ase import Atoms
from ase.build import molecule
from ase.build import bulk
from ase.calculators.lj import LennardJones
from ase.visualize import view
from ase.io import write
from ase.io import read
#from gpaw import GPAW

import numpy as np
import matplotlib.pyplot as plt
#%%
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])

calc = LennardJones()               # create calculator

E = []                              # initialize energy list
rvec = np.linspace(0.8,5.0,100)     # array of distances 100 istället för 2

# loop over distances
for r in rvec:
    co.positions[1][2] = r          # set distance between atoms
    co.set_calculator(calc)         # attach calc to Atoms object
    E0 = co.get_potential_energy()  # calculate energy
    E.append(E0)

plt.figure()
plt.plot(rvec,E)
plt.show()

#%% create atom
h2o = Atoms('H2O', positions=[(0, 0, 0), (0, 0, 1.1), (0, 1.1, 0)])
Na2 = Atoms('Na2', positions=[(0, 0, 0), (0, 0, 1.1)])
view(Na2)

#%% Create molecule
mol = molecule('H2O')
mol.center(vacuum=10.0)

#%% Create chrystal
al = bulk('Al', 'fcc', a=4.05)
print(al.cell)      # automatically the right unitcell

silicon = bulk('Si', 'diamond', a=5.459)
print(silicon.cell)
view(silicon)

#%% Read/write file
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
write('co.xyz',co)                                          # vart skrivs filen!??!
atoms = read('co.xyz') 
view(atoms) 
#%%
#ase build -x diamond Si
#ase build -x diamond Ge
#ase build -x diamond C