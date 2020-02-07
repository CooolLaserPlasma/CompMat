#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TIF320 - Computational materials and molecular physics
Assignment 2
Task 1

"""
from ase import Atom
from ase import Atoms
from ase.build import molecule
from ase.build import bulk
from ase.calculators.lj import LennardJones
from ase.visualize import view
from ase.io import write
from ase.io import read
#from gpaw import GPAW

import numpy as np
import matplotlib.pyplot as plt

from ase import db
#%%
d = 1
cluster = Atoms[('Na', positions=[(0,0,0)]), ('Na', positions=[(0,0,d)])]
view(cluster)