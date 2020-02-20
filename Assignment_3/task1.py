#!/usr/bin/env python3
from gpaw import GPAW
from ase.io.trajectory import Trajectory
from ase.io import read, write
from ase.md.npt import NPT

from ase.units import fs, kB

atoms = read('Na_in_water.xyz')

calc = GPAW(mode='lcao',
            xc='PBE',
            basis='dzp',
            symmetry={'point-group': False},
            charge=1,
            txt='gpawOutput.gpaw-out')

atoms.set_calculator(calc)

dyn = NPT(pfactor=None,
          externalstress=0,
          temperature=350*kB,
          timestep=0.5*fs,
          ttime=20*fs,
          logfile='mdOutput.log')

trajectory = Trajectory('nptDyn.traj', 'w', atoms)
dyn.attach(trajectory.write, interval=1)

#dyn.run(10)
