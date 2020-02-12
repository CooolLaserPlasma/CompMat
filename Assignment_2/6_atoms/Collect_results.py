#!/usr/bin/env python3

from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import write

#########################################################


traj = Trajectory('all_candidates.traj') 	# traj is sorted by energy
atoms = traj[0] 				# Atom configuration with lowest energy

#view(atoms)

# write('least_energy_traj_Na6.xyz', atoms)       # Write atoms object to .xyz file

########################################################
from ase.db import connect
db = connect('gadb.db')

atoms2 = db.select(sort='energy')  # atoms2 is a generator object

atoms2_sorted_list=[]

for atom in atoms2:
        atoms2_sorted_list.append(atom.toatoms())


#print(atoms2_sorted_list)
