#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
from gpaw import GPAW, FermiDirac
from ase.optimize import BFGS
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


def relax(Cluster_to_relax, mod, bas, name):
    # Cluster_to_relax : Atoms object to be optimized
    # mod : GPAW mode parameter e.g. 'lcao', PW(), fd
    # bas : lcao basis e.g. 'dz', 'dzp' etc
    # exc : exchange correlation term
    # name : string with name of txt file and qpw file
    name_of_txt_file = name + '.txt'
    name_of_gpw_file = name + '.gpw'
    
    # Initialize calculator
    calc = GPAW(setups={'Na': '1'}, mode=mod, basis=bas, txt=name_of_txt_file)  	
    Cluster_to_relax.set_calculator(calc)
    dyn = BFGS(Cluster_to_relax)
    dyn.run(fmax=0.02)
    Cluster_to_relax.get_potential_energy()
    #calc.write(name_of_gpw_file)

# OBS lägger bara till occupation, setup, h
# define parameters
mode_of_GPAW = 'lcao'
mode_str = 'lcao'
basis_set = 'dzp'
run_nr = 7
number_of_bands = '-'

# relax structures
relax(Na6_lowest, mode_of_GPAW, basis_set,'struc1_run'+str(run_nr))
relax(Na6_2nd_lowest, mode_of_GPAW, basis_set,'struc2_run'+str(run_nr))


# write result to database
db.write(Na6_lowest, structure=1, GPAWmode=mode_str, run=run_nr, bench='yes-Fermi')
db.write(Na6_2nd_lowest, structure=2, GPAWmode=mode_str, run=run_nr, bench='yes-Fermi')



"""
run 1:
mode='lcao'
xc='PBE'

run 2:
mode=PW(350)
xc='PBE'

run 3:
mode='lcao'
xc='LDA'

run 4:
same parameters as in 'ga.py'
nbands=10
h=0.25
occupations=FermiDirac(0.05)
setups={'Na': '1'}
mode='lcao'
basis='dzp'

run 5:
same as 4 but without Fermi specification

run 6:
same as 4 but without Fermi specification
and without nband

run 7:
same as 4 but without Fermi specification
and without nband
and without h
"""




