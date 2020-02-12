from ase import Atom
from ase import Atoms
from gpaw import GPAW, PW
from ase.visualize import view
#from ase.io.trajectory import Trajectory
#from ase.io import write
from ase.io import read
from ase.db import connect

#%%

#Na6_lowest = read('christmas-tree.xyz')
#Na6_2nd_lowest = read('half-decahedron.xyz')

#%%
#!/usr/bin/env python3


#######################################################

# Create a cluster of Na-atoms in an guess configuration
d = 3

NaCluster = Atoms([Atom('Na', (0,0,0)),
         	   Atom('Na', (d,0,0))])

NaCluster.center(vacuum=5)		# Create surrounding
#view(NaCluster)

#%%
def relax(Cluster_to_relax, mod, name):
    # Cluster_to_relax : Atoms object to be optimized
    # mod : string wit GPAW parameter 'lcao', 'PW'
    # name : string with name of txt file and qpw file
    str1 = name + '.txt'
    str2 = name + '.gpw'
    
    # Initialize calculator
    calc = GPAW(mode=mod, txt=str1)  	
    Cluster_to_relax.set_calculator(calc)
    
    E = Cluster_to_relax.get_potential_energy()
    print(E)
    calc.write(str2)

