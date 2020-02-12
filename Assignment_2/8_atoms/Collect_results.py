#!/usr/bin/env python3

from ase.visualize import view
#from ase.io.trajectory import Trajectory
from ase.io import write
from ase.db import connect

import numpy as np
import matplotlib.pyplot as plt

#%%
db = connect('gadb.db')

# sort db by energy and only include relaxed candidates
atoms = db.select(sort='energy', relaxed=1)         # atoms is a generator object

atoms_sorted_list=[] # list that will contain atoms objects sorted by energy 

for atom in atoms:
        atoms_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms_sorted_list))

for i in range(len(atoms_sorted_list)):
    energy[i] = atoms_sorted_list[i].get_potential_energy()

diffs = abs(energy - np.roll(energy, -1)) # diffrences in energy between adjecent configurations
diffs = diffs[0:-1]  # remove last entry since this is the difference between lowest and highest energy

#%%
# plot energy differences
"""
where the energy difference is larger than about 0.03 there seems to be a new structure
"""
fig, ax = plt.subplots(1, 1)
ax.plot(diffs)
ax.set_xlabel('"Canditade index"', fontsize=15)
ax.set_ylabel('Energy difference (eV)', fontsize=15)
ax.plot([0, 100], [0.03, 0.03], '--')
#fig.savefig('Energydifferences_Na7.png')


#%%

diffs_list = [i for i in diffs]   # make diffs as a list 
index = [diffs_list.index(i) for i in diffs_list if i > 0.03] # pick out indices for diffs > 0.03
"""
index contains indices of "atoms_sorted_list" where a new configuration is
index[0]+1 will get the first new configuration
"""
lowest_energy_structure = atoms_sorted_list[0]
E_lowest = lowest_energy_structure.get_potential_energy()

snd_lowest_energy_structure = atoms_sorted_list[index[0]+1]
E_snd_lowest = snd_lowest_energy_structure.get_potential_energy()

print('Energy for most stabel structure: ' + str(E_lowest))
print('Energy for 2nd most stabel structure: ' + str(E_snd_lowest))

#%%
#view(lowest_energy_structure)
#view(snd_lowest_energy_structure)

# save structures as .xyz files
#write('lowest_energy_struc_Na7.xyz', lowest_energy_structure)
#write('2nd_lowest_energy_struc_Na7.xyz', snd_lowest_energy_structure)


