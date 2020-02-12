#!/usr/bin/env python3

from ase import Atom
from ase import Atoms
#from ase.build import molecule
#from ase.build import bulk
#from ase.calculators.lj import LennardJones
from ase.visualize import view
#from ase.io import write
#from ase.io import read
from gpaw import GPAW

#import numpy as np
#import matplotlib.pyplot as plt

#from ase import db
############################################	
#Na = Atom('Na')


# Create a cluster of Na-atoms in an guess configuration
d = 2
NaCluster = Atoms([Atom('Na', (0,0,0)),
		   Atom('Na', (d,0,0)),
		   Atom('Na', (0,d,0)),
		   Atom('Na', (0,0,d)),
		   Atom('Na', (d,d,0)),
		   Atom('Na', (0,d,d))])

NaCluster.center(vacuum=10)
view(NaCluster)

# initialize a calculator


############################################
# Code from GPAW-manual
#d = 0.74
#a = 6.0

#atoms = Atoms('H2',
#              positions=[(0, 0, 0),
#                         (0, 0, d)],
#              cell=(a, a, a))
#atoms.center()

#calc = GPAW(nbands=2, txt='h2.txt')
#atoms.set_calculator(calc)
#print(atoms.get_forces())

# trying without specifying a cell and not "atoms.center()"
# DOES NOT WORK
# same bu with "atoms.center()" gives
# ValueError: GPAW requires 3 lattice vectors.  Your system has 0.
# same but with "atom.center(vacuum=5)" seems to work!! :)
#d = 0.74
#a = 6.0

#atoms = Atoms('H2',
#              positions=[(0, 0, 0),
#                         (0, 0, d)])
#atoms.center(vacuum=5)

#calc = GPAW(nbands=2, txt='h2.txt')
#atoms.set_calculator(calc)

#E = atoms.get_potential_energy()
#print(E)
#calc.write('h2.gpw')   # added this from the code  below

###################################################
# part of code from GPAW documantation: "Calculation of atomization energies"

#a = 10.  # Size of unit cell (Angstrom)
#c = a / 2

# Hydrogen atom:
#atom = Atoms('H',
#             positions=[(c, c, c)],
#             magmoms=[0],
#             cell=(a, a + 0.0001, a + 0.0002))  # Break cell symmetry

# gpaw calculator:
#calc = GPAW(mode=PW(),
#            xc='PBE',
#            hund=True,
#            eigensolver='rmm-diis',  # This solver can parallelize over bands
#            occupations=FermiDirac(0.0, fixmagmom=True),
#            txt='H.out',
#            )
#atom.set_calculator(calc)

#e1 = atom.get_potential_energy()
#calc.write('H.gpw')

########################################################

	
