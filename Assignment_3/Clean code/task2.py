#!/usr/bin/env python3
from ase.io.trajectory import Trajectory

import numpy as np
import matplotlib.pyplot as plt

# load files
traj_7ps = Trajectory('NaCluster24.traj')
traj_2ps = Trajectory('nptDyn.traj')  

log_7ps = open('NaCluster24.log', 'r')
t_7ps, E_tot_7ps, E_pot_7ps, E_kin_7ps, T_7ps = np.loadtxt(log_7ps, 
                                                            skiprows=1, 
                                                            unpack=True)

log_2ps = open('mdOutput.log', 'r')
t_2ps, E_tot_2ps, E_pot_2ps, E_kin_2ps, T_2ps = np.loadtxt(log_2ps, 
                                                            skiprows=1, 
                                                            unpack=True)

# define some functions
def get_r(trajectory, equilibrium_frame):
    # gets distances between the Na atom and every oxygen atom for all
    # equilibrated snapshots (frames)
    # trajectory : trajectory to extract distances from
    # eq_frame   : index of equilibrated frame
    
    indices = np.arange(0,24)       # indices of all oxygen atoms (0,..,23)
                                    # the Na atom has the last index (-1)
    r_vect = []
    nframes = len(trajectory) - equilibrium_frame 

    for i in range(nframes):        
        atoms = trajectory[equilibrium_frame + i]
        r_vect.append(atoms.get_distances(-1, indices, mic=True))

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

# set index of first equilibrated frame
eq_7ps = 1000
eq_2ps = 1000

# get the distances between the Na atom and every oxygen atom 
r_7ps = get_r(traj_7ps, eq_7ps)    
r_2ps = get_r(traj_2ps, eq_2ps)

# get radial distribution functions
number_of_bins = 1000  
[g_7ps, R_7ps] = RDF(r_7ps, number_of_bins, traj_7ps, eq_7ps)
[g_2ps, R_2ps] = RDF(r_2ps, number_of_bins, traj_2ps, eq_2ps)


# plot radial distribution functions
fig, ax = plt.subplots(1, 2, figsize=(12, 4))

ind_7ps = 225 # index of first minima
ax[0].plot(R_7ps[0:-1], g_7ps)
ax[0].set_xlabel('r [Å]', fontsize=20)
ax[0].set_ylabel('g(r)', fontsize=20)
ax[0].plot([R_7ps[ind_7ps], R_7ps[ind_7ps]], [0, 5], '--')
ax[0].set_title('7ps', fontsize=20)
ax[0].tick_params(axis='both', which='major', labelsize=15)

ind_2ps = 201
ax[1].plot(R_2ps[0:-1], g_2ps)
ax[1].set_xlabel('r [Å]', fontsize=20)
ax[1].set_ylabel('g(r)', fontsize=20)
ax[1].plot([R_2ps[ind_2ps], R_2ps[ind_2ps]], [0, 6], '--')
ax[1].set_title('2ps', fontsize=20)
ax[1].tick_params(axis='both', which='major', labelsize=15)

plt.tight_layout()

# Calculate coordination numbers
coordnr_7ps = coordination_nr(g_7ps, R_7ps, ind_7ps, traj_7ps) 
coordnr_2ps = coordination_nr(g_2ps, R_2ps, ind_2ps, traj_2ps) 

print('Coordination numbers')
print('7ps: ' + str(coordnr_7ps))
print('2ps: ' + str(coordnr_2ps))