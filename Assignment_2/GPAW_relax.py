#!/usr/bin/env python3

from ase import Atom
from ase import Atoms
from gpaw import GPAW, PW

#######################################################

# Create a cluster of Na-atoms in an guess configuration
d = 2
#NaCluster = Atoms([Atom('Na', (0,0,0)),
#         	   Atom('Na', (d,0,0)),
#            	   Atom('Na', (0,d,0)),
#            	   Atom('Na', (0,0,d)),
#            	   Atom('Na', (d,d,0)),
#            	   Atom('Na', (0,d,d))])

NaCluster = Atoms([Atom('Na', (0,0,0)),
                  Atom('Na', (d,0,0))])

NaCluster.center(vacuum=5)		# Create surrounding

# Initialize calculator
calc = GPAW(mode=PW(), txt='Na6.txt')  	# complaining on this, try import PW
NaCluster.set_calculator(calc)

E = NaCluster.get_potential_energy()
print(E)
calc.write('Na6.gpw') 
