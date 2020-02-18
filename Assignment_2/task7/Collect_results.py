#!/usr/bin/env python3
from ase.visualize import view
from ase.io import write
from ase.db import connect

import numpy as np
import matplotlib.pyplot as plt

# Na6
db = connect('gadb_Na6.db')

# sort db by energy and only include relaxed candidates
atoms = db.select(sort='energy', relaxed=1)

atoms_sorted_list=[] # list that will contain atoms objects sorted by energy 

for atom in atoms:
        atoms_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms_sorted_list))

for i in range(len(atoms_sorted_list)):
    energy[i] = atoms_sorted_list[i].get_potential_energy()

diffs = abs(energy - np.roll(energy, -1)) # diffrences in energy between adjecent configurations
diffs = diffs[0:-1]  # remove last entry since this is the difference between lowest and highest energy


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

E_diff = abs(E_lowest-E_snd_lowest)

print('Energy difference for Na6: ' + str(E_diff))

# Na7
db = connect('gadb_Na7.db')

# sort db by energy and only include relaxed candidates
atoms = db.select(sort='energy', relaxed=1)

atoms_sorted_list=[] # list that will contain atoms objects sorted by energy 

for atom in atoms:
        atoms_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms_sorted_list))

for i in range(len(atoms_sorted_list)):
    energy[i] = atoms_sorted_list[i].get_potential_energy()

diffs = abs(energy - np.roll(energy, -1)) # diffrences in energy between adjecent configurations
diffs = diffs[0:-1]  # remove last entry since this is the difference between lowest and highest energy


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

E_diff = abs(E_lowest-E_snd_lowest)

print('Energy difference for Na7: ' + str(E_diff))

# Na8 
db = connect('gadb_Na8.db')

# sort db by energy and only include relaxed candidates
atoms = db.select(sort='energy', relaxed=1)

atoms_sorted_list=[] # list that will contain atoms objects sorted by energy 

for atom in atoms:
        atoms_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms_sorted_list))

for i in range(len(atoms_sorted_list)):
    energy[i] = atoms_sorted_list[i].get_potential_energy()

diffs = abs(energy - np.roll(energy, -1)) # diffrences in energy between adjecent configurations
diffs = diffs[0:-1]  # remove last entry since this is the difference between lowest and highest energy


diffs_list = [i for i in diffs]   # make diffs as a list 
index = [diffs_list.index(i) for i in diffs_list if i > 0.03] # pick out indices for diffs > 0.03
"""
index contains indices of "atoms_sorted_list" where a new configuration is
index[0]+1 will get the first new configuration
"""

lowest_energy_structure = atoms_sorted_list[0]
E_lowest = lowest_energy_structure.get_potential_energy()

snd_lowest_energy_structure = atoms_sorted_list[index[1]+1] 
E_snd_lowest = snd_lowest_energy_structure.get_potential_energy()

E_diff = abs(E_lowest-E_snd_lowest)

print('Energy difference for Na8: ' + str(E_diff))
