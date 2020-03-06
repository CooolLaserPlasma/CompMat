#!/usr/bin/env python3
from ase.io.trajectory import Trajectory

import numpy as np
import matplotlib.pyplot as plt

# load files
traj = Trajectory('cluster24.traj')  

log = open('cluster24.log', 'r')
t, E_tot, E_pot, E_kin, T = np.loadtxt(log, 
                                       skiprows=1, 
                                       unpack=True)


# define some functions
def get_r(trajectory, equilibrium_frame):
    # gets distances between the Na atom and every oxygen atom for all
    # equilibrated snapshots (frames)
    # trajectory : trajectory to extract distances from
    # eq_frame   : index of equilibrated frame
    
    indices = np.arange(0,23)       # indices of all oxygen atoms
                                    # except the first one 
    r_vect = []
    nframes = len(trajectory) - equilibrium_frame 

    for i in range(nframes):        
        atoms = trajectory[equilibrium_frame + i]
        r_vect.append(atoms.get_distances(23, indices, mic=True))

    r_array = np.asarray(r_vect)
    r_array = r_array.reshape(1,len(indices)*nframes)
    
    return r_array

def RDF(radial_distances, nbins, trajectory, equilibrium_frame):
    # Gives the radial distribution function from a vector whit all relevant
    # distances
    # radial_distances : vector with radial distances
    # nbins            : number of bins for the histogram 

    nframes = len(trajectory) - equilibrium_frame
    volume = trajectory[-1].get_volume() 
    rho =  24/volume                             

    [antal, R] = np.histogram(radial_distances, bins=nbins) 
    dr = R[1]-R[0]
    radial_distribution_function = antal/(nframes*dr*rho*4*np.pi*R[0:-1]**2)
    return [radial_distribution_function, R]

def coordination_nr(radial_dist_func, distance_vector, first_minimum, trajectory):
    dr = distance_vector[1]-distance_vector[0]
    volume = trajectory[-1].get_volume()
    rho = 24/volume
    coordination_number = sum(4*np.pi*distance_vector[0:first_minimum]**2*radial_dist_func[0:first_minimum]*dr*rho)
    return coordination_number

# set index of first equilibrium frame
equilibrium = 1000
# get the distances between the Na atom and every oxygen atom  
r = get_r(traj, equilibrium)    

# get radial distribution functions 
[g, R] = RDF(r, 500, traj, equilibrium)

# plot radial distribution functions
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.plot(R[0:-1], g)
ax.set_xlabel('r [Å]', fontsize=20)
ax.set_ylabel('g(r)', fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

plt.tight_layout()
