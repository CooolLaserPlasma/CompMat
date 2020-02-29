#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
#from gpaw import GPAW, PW
from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import read, write

import numpy as np
import matplotlib.pyplot as plt
#%%
traj = Trajectory('NaCluster24.traj')  # 7ps AIMD simulation of 24 water molecules + 1 Na ion


log = open('NaCluster24.log', 'r')
t, E_tot, E_pot, E_kin, T = np.loadtxt(log, skiprows=1, unpack=True)


mean_T = sum(T)/len(T)

#%%
# check for thermal equilibrium 
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(t, T)
ax[0].plot([0, 6],[mean_T, mean_T])
ax[0].set_xlabel('time [ps]', fontsize=15)
ax[0].set_ylabel('temperature [K]', fontsize=15)

t_max = 3000
ax[1].plot(t[0:t_max], T[0:t_max])
ax[1].plot([0, 1.5],[mean_T, mean_T])
ax[1].set_xlabel('time [ps]', fontsize=15)
ax[1].set_ylabel('temperature [K]', fontsize=15)
#fig.savefig('temp_vs_time.png')

timestep = t[1]-t[0]   # timestep used in the AIMD simulation
n = int(1/timestep)

#%%
# get distances etc
# OBS! dont forget periodic bcs:)

# for each snashot (do only include equillibrated timesteps!)
# how do I get atoms from trajectory?

atoms = traj[1000]
#symbols = atoms.get_chemical_symbols()
#atom = atoms[-1]                        # the Na atom has the last index
indices = np.arange(0,24)               # indices of all O atoms
r_vect = atoms.get_distances(-1, indices, mic=True) # distance from Na atom to each O atom

#%%
r_vect = []
# range to specify relevant trajectories
nframes = len(traj)

for i in range(nframes):
    atoms = traj[i] 
    r_vect.append(atoms.get_distances(-1, indices, mic=True))

#%%
r_array = np.asarray(r_vect)
r_array.reshape(1,24*nframes)

volume = atoms.get_volume() 
rho =  24/volume                             # should include Naq ion aswell!?

[antal, R] = np.histogram(r_array, bins=100) # antal per spherical shell at distance r
dr = R[1]-R[0]
g = antal/(nframes*dr*rho*4*np.pi*R[0:-1]**2)

fig, ax = plt.subplots(1, 1)
ax.plot(R[0:-1], g)
ax.set_xlabel('r', fontsize=15)
ax.set_ylabel('g(r)', fontsize=15)

#fig.savefig('Energydifferences.png')
#%%
# calculate coordination number
firstmin = 20           # OBS 20 är för bins=100
#4 * pi * x^2 * RDF(x) and multiply with the average density rho = len(atoms) / atoms.get_volume()
coordnr = sum(4*np.pi*R[0:firstmin]**2*g[0:firstmin]*dr*rho)
print(coordnr)

