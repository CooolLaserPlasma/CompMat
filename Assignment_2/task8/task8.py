#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
from gpaw import GPAW, PW
from ase.visualize import view
#from ase.io.trajectory import Trajectory
#from ase.io import write
from ase.io import read
from ase.db import connect


# Read structures to be investigated
Na6_lowest = read('christmas-tree.xyz')
Na6_2nd_lowest = read('half-decahedron.xyz')

# Connect to database
db = connect('t8db.db')


def relax(Cluster_to_relax, mod, name):
    # Cluster_to_relax : Atoms object to be optimized
    # mod : GPAW mode parameter e.g. 'lcao', PW()
    # name : string with name of txt file and qpw file
    name_of_txt_file = name + '.txt'
    name_of_gpw_file = name + '.gpw'
    
    # Initialize calculator
    calc = GPAW(mode=mod, txt=name_of_txt_file)  	
    Cluster_to_relax.set_calculator(calc)
    
    #E = Cluster_to_relax.get_potential_energy()
    calc.write(name_of_gpw_file)


# define parameters
mode_of_GPAW = 'lcao'
run_nr = 1

# relax structures
relax(Na6_lowest, mode_of_GPAW, 'struc1_run'+str(run_nr))
relax(Na6_2nd_lowest, mode_of_GPAW, 'struc2_run'+str(run_nr))


# write result to database
db.write(Na6_lowest, structure=1, GPAWmode=mode_of_GPAW, run=run_nr)
db.write(Na6_2nd_lowest, structure=2, GPAWmode=mode_of_GPAW, run=run_nr)


"""
run 1:
mode = 'lcao'



"""




