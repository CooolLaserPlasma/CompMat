#!/usr/bin/env python3

from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import write

#########################################################


traj = Trajectory('all_candidates.traj') 	# traj is sorted by energy
atoms = traj[0] 				# Atom configuration with lowest energy

view(atoms)

# write('least_energy_traj_Na6.xyz', atoms)       # Write atoms object to .xyz file
