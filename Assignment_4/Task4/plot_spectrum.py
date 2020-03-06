#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import read, write

import numpy as np
import matplotlib.pyplot as plt
#%%
# load files
spectrum_x = open('be_spectrum_x.dat', 'r')
spectrum_y = open('be_spectrum_y.dat', 'r')
spectrum_z = open('be_spectrum_z.dat', 'r')

om, S_x, y, z = np.loadtxt(spectrum_x, 
                           skiprows=6, 
                           unpack=True)

om, x, S_y, z = np.loadtxt(spectrum_y, 
                           skiprows=6, 
                           unpack=True)

om, x, y, S_z = np.loadtxt(spectrum_z, 
                           skiprows=6, 
                           unpack=True)

#%%
start = 0
stop = 120
# check for thermal equilibrium 
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.plot(om[start:stop], S_x[start:stop], label='x')
ax.plot(om[start:stop], S_y[start:stop], label='y')
ax.plot(om[start:stop], S_z[start:stop], label='z')
ax.set_xlabel('Energy [eV]', fontsize=20)
ax.set_ylabel('Photoabsorption', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(fontsize=20)

plt.tight_layout()
#fig.savefig('photoabsorption.png')

#%%
mean_S = (S_x + S_y + S_z)/3

start = 0
stop = 120
# check for thermal equilibrium 
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.plot(om[start:stop], mean_S[start:stop])
ax.set_xlabel('Energy [eV]', fontsize=20)
ax.set_ylabel('mean photoabsorption', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

plt.tight_layout()
fig.savefig('mean_photoabsorption.png')

#%%
# plot dipole moments
dm_x = open('be_dm_x.dat', 'r')
dm_y = open('be_dm_y.dat', 'r')
dm_z = open('be_dm_z.dat', 'r')

t, norm_x, x, y, z = np.loadtxt(dm_x, 
                                skiprows=2, 
                                unpack=True)
dmx_norm = np.sqrt(x + y + z)

t, norm_y, x, y, z = np.loadtxt(dm_y, 
                                skiprows=2, 
                                unpack=True)
dmy_norm = np.sqrt(x + y + z)

t, norm_z, x, y, z = np.loadtxt(dm_z, 
                                skiprows=2, 
                                unpack=True)
dmz_norm = np.sqrt(x + y + z)
#%%
norm = norm_x + norm_y + norm_z
start = 0
stop = len(t)
# check for thermal equilibrium 
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.plot(t[start:stop], dmx_norm[start:stop], label='x')
ax.plot(t[start:stop], dmy_norm[start:stop], label='y')
ax.plot(t[start:stop], dmz_norm[start:stop], label='z')
ax.set_xlabel('Energy [eV]', fontsize=20)
ax.set_ylabel('Photoabsorption', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.legend(fontsize=20)

plt.tight_layout()
#fig.savefig('photoabsorption.png')