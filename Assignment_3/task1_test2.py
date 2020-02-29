#!/usr/bin/env python3
from ase import Atom
from ase import Atoms
#from gpaw import GPAW, PW
from ase.visualize import view
from ase.io.trajectory import Trajectory
from ase.io import read, write

import numpy as np
import matplotlib.pyplot as plt
#%%
# load files
traj_7ps = Trajectory('NaCluster24.traj')  # 7ps AIMD simulation of 24 water molecules + 1 Na ion

log_7ps = open('NaCluster24.log', 'r')
t_7ps, E_tot_7_ps, E_pot_7ps, E_kin_7ps, T_7ps = np.loadtxt(log_7ps, 
                                                            skiprows=1, 
                                                            unpack=True)

#%%
# check for thermal equilibrium 
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(t_7ps, T_7ps)
ax[0].set_xlabel('time [ps]', fontsize=15)
ax[0].set_ylabel('temperature [K]', fontsize=15)

t_max = 3000
ax[1].plot(t_7ps[0:t_max], T_7ps[0:t_max])
ax[1].set_xlabel('time [ps]', fontsize=15)
ax[1].set_ylabel('temperature [K]', fontsize=15)
#fig.savefig('temp_vs_time.png')

#%% 
# define some functions
def get_r(trajectory, equilibrium_frame):
    # gets distances between the Na atom and every oxygen atom for all
    # equilibrated snapshots (frames)
    # trajectory : trajectory to extract distances from
    # eq_frame   : index of equilibrated frame
    
    indices = np.arange(0,24)       # indices of all oxygen atoms
                                    # the Na atom has the last index (-1)
    r_vect = []
    nframes = len(trajectory) - equilibrium_frame 

    for i in range(nframes):        
        atoms = trajectory[equilibrium_frame + i] 
        r_vect.append(atoms.get_distances(-1, indices, mic=True))

    r_array = np.asarray(r_vect)
    r_array = r_array.reshape(1,24*nframes)
    
    return r_array

def RDF(radial_distances, nbins, trajectory, equilibrium_frame):
    # Gives the radial distribution function from a vector whit all relevant
    # distances
    # radial_distances : vector with radial distances
    # nbins            : number of bins for the histogram 

    nframes = len(trajectory) - equilibrium_frame
    volume = traj_7ps[-1].get_volume() 
    rho =  24/volume                             # should include Naq ion aswell!?

    [antal, R] = np.histogram(radial_distances, bins=nbins) # antal per spherical shell at distance r
    dr = R[1]-R[0]
    radial_distribution_function = antal/(nframes*dr*rho*4*np.pi*R[0:-1]**2)
    return [radial_distribution_function, R]

def coordination_nr(radial_dist_func, distance_vector, first_minimum, trajectory):
    dr = distance_vector[1]-distance_vector[0]
    volume = trajectory[-1].get_volume()
    rho = 24/volume
    coordination_number = sum(4*np.pi*distance_vector[0:first_minimum]**2*radial_dist_func[0:first_minimum]*dr*rho)
    return coordination_number
#%% (This step takes a while...)
# get the distances between the Na atom and every oxygen atom
    
r_7ps = get_r(traj_7ps, 13500)    


#%%
# get radial distribution functions 
[g_7ps, R_7ps] = RDF(r_7ps, 100, traj_7ps, 13500)



# plot radial distribution functions
fig, ax = plt.subplots(1, 1)
ax.plot(R_7ps[0:-1], g_7ps)
ax.set_xlabel('r', fontsize=15)
ax.set_ylabel('g(r)', fontsize=15)

#fig.savefig('Energydifferences.png')
#%%
# Calculate coordination numbers
    
coordnr_7ps = coordination_nr(g_7ps, R_7ps, 20, traj_7ps) # 20 for 100 bins
print(coordnr_7ps)
