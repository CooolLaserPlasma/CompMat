#!/usr/bin/env python3

from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import write
import numpy as np
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

atoms2_sorted_list=[]

for atom in atoms2:
        atoms2_sorted_list.append(atom.toatoms())


# get diffrences in energy)
energy = np.zeros(len(atoms2_sorted_list))

for i in range(len(atoms2_sorted_list)):
    energy[i] = atoms2_sorted_list[i].get_potential_energy()


print(len(energy))
diffs = energy - np.roll(energy, -1)
#print(diffs)



Na6_lowest = atoms2_sorted_list[2]
#view(Na6_lowest)



###############
Na6_lowest2 = db.get('id=74').toatoms()


#view(Na6_lowest2)
