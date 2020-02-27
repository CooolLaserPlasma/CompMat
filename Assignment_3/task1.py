#!/usr/bin/env python3
from gpaw import GPAW
from ase.io.trajectory import Trajectory
from ase.io import read, write
from ase.md.npt import NPT
from ase.visualize import view
from ase.units import fs, kB
import numpy as np

atoms = read('Na_in_water.xyz')
#view(atoms)

#com = atoms.get_momenta().sum(0)
#print(np.sqrt(np.sum(com*com)))

calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            symmetry={'point_group': False},
            charge=1,
            txt='gpawOutput.gpaw-out')

atoms.set_calculator(calc)

dyn = NPT(atoms,
	  pfactor=None,
          externalstress=0,
          temperature=350*kB,
          timestep=0.5*fs,
          ttime=20*fs,
          logfile='mdOutput.log')

trajectory = Trajectory('nptDyn.traj', 'w', atoms)
dyn.attach(trajectory.write, interval=1)


dyn.run(4000)
