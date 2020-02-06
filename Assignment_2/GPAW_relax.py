#!/usr/bin/env python3

from ase import Atom
from ase import Atoms
from gpaw import GPAW, PW
from ase.visualize import view

#######################################################

# Create a cluster of Na-atoms in an guess configuration
d = 3

NaCluster = Atoms([Atom('Na', (0,0,0)),
         	   Atom('Na', (d,0,0)),
            	   Atom('Na', (0,d,0)),
            	   Atom('Na', (0,0,d)),
            	   Atom('Na', (d,d,0)),
            	   Atom('Na', (0,d,d))])

NaCluster.center(vacuum=5)		# Create surrounding
#view(NaCluster)


# Initialize calculator
calc = GPAW(mode='lcao', txt='Na6.txt')  	
NaCluster.set_calculator(calc)

E = NaCluster.get_potential_energy()
print(E)
calc.write('Na6.gpw')



"""
1st run on hebbe:
d = 2
vacuum = 5
mode = PW()
did not converge

2nd run on hebbe:
d = 3
vacuum = 5
mode = PW()

3rd run on hebbe:
d = 3
vacuum = 5
mode = 'lcao'

""" 
