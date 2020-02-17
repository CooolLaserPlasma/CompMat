#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
from ase.io import read


# Read structures
Na6_lowest = read('christmas-tree.xyz')
Na6_2nd_lowest = read('half-decahedron.xyz')

# Get potential energies
E_christmas = Na6_lowest.get_potential_energy()
E_decahedron = Na6_2nd_lowest.get_potential_energy()

print('Potential energy of the two most stable structures before relaxation')
print('Christmas-tree: '+str(E_christmas))
print('Half_decahedron: '+str(E_decahedron))

