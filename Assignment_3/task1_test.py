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
traj = Trajectory('cluster24.traj')  # 6 ps AIMD simulation of 24 water molecules


log = open('cluster24.log', 'r')
t, E_tot, E_pot, E_kin, T = np.loadtxt(log, skiprows=1, unpack=True)


mean_T = sum(T)/len(T)

#%%
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
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(t, E_tot)
ax[0].set_xlabel('time [ps]', fontsize=15)
ax[0].set_ylabel('energy [eV]', fontsize=15)

t_max = 3000
ax[1].plot(t[0:t_max], E_tot[0:t_max])
ax[1].set_xlabel('time [ps]', fontsize=15)
ax[1].set_ylabel('energy [eV]', fontsize=15)

#%%

config_nr = n               # n+7 has temp 347.5 (T_mean = 345.1)

atoms = traj[config_nr]
print(T[config_nr])
#view(atoms)

com = atoms.get_momenta().sum(0)
print(np.sqrt(np.sum(com*com)))

#atoms_Na = read('water_Na.xyz')
#view(atoms_Na)

#%%
    
#Na_ion = Atoms(Atom('Na', charge=1), (0,0,0))
#write('Na_ion.xyz', Na_ion)



