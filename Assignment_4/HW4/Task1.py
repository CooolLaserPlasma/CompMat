from gpaw import GPAW
from ase.io import read
from gpaw.lrtddft import LrTDDFT
from gpaw.lrtddft import photoabsorption_spectrum
import numpy as np
import matplotlib.pyplot as plt


def calc_omega(gpw_file):
	calc = GPAW(gpw_file)	
	dE = 6  # maximal Kohn-Sham transition energy to consider in eV
	lr = LrTDDFT(calc, xc='LDA', energy_range=dE)
	lr.write('lr_dE.dat.gz')
	lr.diagonalize()
	# write the spectrum to the data file
	photoabsorption_spectrum(lr, 'spectrum_w,06eV.dat', # data file name
                         	 width=0.06)                # width in eV

def init_bands():
	atoms = read('Na-tddft/Na8.xyz')
	atoms.center(vacuum=8)

	calc = GPAW(mode = 'fd',
		    xc = 'LDA',
		    setups = {'Na': '1'}, # Use a single - valence electron setup
		    h = 0.3,
		    nbands = 0, # Include only occupied states
	            txt = 'f.gpaw-out')

	atoms.set_calculator(calc)
	atoms.get_potential_energy() # Performs the calculation on the system !
	calc.write('0bands.gpw', 'all') # Writes everything it can to file !

def band_calc(bands):
	calc = GPAW('0bands.gpw')
	calc.set(nbands = bands,
		 fixdensity = True, # Fixes density in any new calculations
		 convergence={'bands': -10}) # converge all but the last 10
	atoms = calc.get_atoms()
	atoms.get_potential_energy() # Recalculate with new parameters !
	calc.write('bands'+str(bands)+'.gpw', 'all') # Writes everything it can to file !

def print_eig(gpw_filename):
	calc = GPAW(gpw_filename,txt=None)
	#atoms = calc.get_atoms()
	eig = calc.get_eigenvalues()
	print('eigenvalues:\n',eig)
	print('energy difference highest occupied and highest unoccupied: ' 
	      + str(eig[-11]-eig[3]))

def plot_spectrum(filename):
	data  = np.genfromtxt(filename)
	plt.figure()
	plt.plot(data[:,1])
	plt.show()
	print(data)

if __name__=='__main__':
	#band_calc(110)
	#print_eig('bands110.gpw')
	#calc_omega('bands110.gpw')
	plot_spectrum('spectrum_w,06eV.dat')
