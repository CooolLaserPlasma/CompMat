#!/usr/bin/env python3

#from ase import Atom
#from ase import Atoms
#from gpaw import GPAW, PW
from ase.visualize import view

#from ase.db import connect
from ase.io.trajectory import Trajectory
from ase.io import write
#########################################################

#db = connect('gadb.db')
#E = db.select(0 ,sort='enegry')
#print(E)



traj = Trajectory('all_candidates.traj') 	# traj is sorted by energy
atoms = traj[0] 				# Atom configuration with lowest energy

#E = atoms.get_potential_energy()
#print(E)
view(atoms)

write('least_energy_traj_Na6.xyz', atoms)
