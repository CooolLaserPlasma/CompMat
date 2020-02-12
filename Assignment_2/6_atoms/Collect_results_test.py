#!/usr/bin/env python3

from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import write
import numpy as np

import matplotlib.pyplot as plt
#%%
#########################################################


#traj = Trajectory('all_candidates.traj') 	# traj is sorted by energy
#atoms = traj[0] 				# Atom configuration with lowest energy

#view(atoms)

# write('least_energy_traj_Na6.xyz', atoms)          # Write atoms object to .xyz file

########################################################
from ase.db import connect
db = connect('gadb.db')

# sort by energy and only include relaxed candidates
atoms2 = db.select(sort='energy', relaxed=1)         # atoms2 is a generator object

atoms2_sorted_list=[] # list that will contain atoms objects sorted by energy 

for atom in atoms2:
        atoms2_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms2_sorted_list))

for i in range(len(atoms2_sorted_list)):
    energy[i] = atoms2_sorted_list[i].get_potential_energy()


diffs = abs(energy - np.roll(energy, -1))
diffs = diffs[0:-1]  # remove last entr since this is the difference between lowest and highest energy

#%%
"""
where the energy difference is larger than about 0.03 there seems to be a new structure
"""
fig, ax = plt.subplots(1, 1)
ax.plot(diffs)
ax.set_xlabel('"Canditade index"', fontsize=15)
ax.set_ylabel('Energy difference (eV)', fontsize=15)
ax.plot([0, 100], [0.03, 0.03], '--')
#fig.savefig('Energydifferences.png')


#%%
# try with lists instead??
# redo diffs as a list: diffs2 = [i for i in diffs]
# pick out indices: index = [diffs2.index(i) for i in diffs2 if i > 0.02]

diffs_list = [i for i in diffs]
index = [diffs_list.index(i) for i in diffs_list if i > 0.03] 
"""
index contains indices of "atoms2_sorted_list" where a new configuration is
index[0]+1 will get the first new configuration
"""
Config1 = atoms2_sorted_list[index[1]]
Config2 = atoms2_sorted_list[index[1]+1]


#%%
Na6_lowest = atoms2_sorted_list[2]
#view(Na6_lowest)



###############
Na6_lowest2 = db.get('id=74').toatoms()


#view(Na6_lowest2)


