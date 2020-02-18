#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
from gpaw import GPAW, FermiDirac, PW
from ase.optimize import BFGS
from ase.visualize import view
#from ase.io.trajectory import Trajectory
#from ase.io import write
from ase.io import read
from ase.db import connect
#from gpaw import setup_paths
#setup_paths.insert(0, '/cephyr/users/hallborn/Hebbe/HW2/CompMat/Assignment_2/task8')

# Read structures to be investigated
Na6_lowest = read('christmas-tree.xyz')
Na6_2nd_lowest = read('half-decahedron.xyz')


# Connect to database
db = connect('t8db.db')


def relax(Cluster_to_relax, mod, bas, exc, name):
    # Cluster_to_relax : Atoms object to be optimized
    # mod : GPAW mode parameter e.g. 'lcao', PW(), fd
    # bas : lcao basis e.g. 'dz', 'dzp' etc
    # exc : exchange correlation term
    # name : string with name of txt file and qpw file
    name_of_txt_file = name + '.txt'
    name_of_gpw_file = name + '.gpw'
    
    # Initialize calculator
    calc = GPAW(setups={'Na': '1'}, h=0.25, occupations=FermiDirac(0.05), mode=mod, basis=bas, xc=exc, txt=name_of_txt_file)  	
    Cluster_to_relax.set_calculator(calc)
    dyn = BFGS(Cluster_to_relax)
    dyn.run(fmax=0.02)
    Cluster_to_relax.get_potential_energy()
    #calc.write(name_of_gpw_file)


# define parameters
mode_of_GPAW = PW(200)
mode_str = 'PW(200)'
basis_set = 'dzp'
exchange = 'PBE'
fermi_dirac = 0.05
run_nr = 8

# relax structures
relax(Na6_lowest, mode_of_GPAW, basis_set, exchange, 'struc1_run'+str(run_nr))
relax(Na6_2nd_lowest, mode_of_GPAW, basis_set, exchange, 'struc2_run'+str(run_nr))


# write result to database
db.write(Na6_lowest, structure=1, GPAWmode=mode_str, basis=basis_set, xc=exchange, fermi=fermi_dirac, run=run_nr)
db.write(Na6_2nd_lowest, structure=2, GPAWmode=mode_str, basis=basis_set, xc=exchange, fermi=fermi_dirac, run=run_nr)



"""
run 1:
setups={'Na': '1'}
h=0.25
mode='lcao'
basis='dzp'
xc='PBE'

run 2:
setups={'Na': '1'}
h=0.25
mode='lcao'
basis='sz(dzp)'
xc='PBE'

run 3:
setups={'Na': 'paw'}
h=0.25
mode='lcao'
basis='dzp'
xc='LDA'

run 4:
setups={'Na': 'paw'}
h=0.25
mode='lcao'
basis='sz(dzp)'
xc='LDA'

run 5:
setups={'Na': 'paw'}
h=0.25
mode='lcao'
basis='qztp'
xc='PBE'

run 6:
setups={'Na': '1'}
h=0.25
mode=PW(350)
basis='dzp'
xc='PBE'

run 7:
setups={'Na': '1'}
h=0.25
mode=PW(500)
basis='dzp'
xc='PBE'

run 8:
setups={'Na': '1'}
h=0.25
mode=PW(200)
basis='dzp'
xc='PBE'

run 8:  (yes double 8)
setups={'Na': '1'}
h=0.25
mode=PW(200)
basis='dzp'
xc='PBE'
occupations=FermiDirac(0.05)
"""




