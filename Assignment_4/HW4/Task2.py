import sys
sys.path.append('./Na-tddft/')
from helper import dump_data, fold
import numpy as np
from gpaw.lrtddft import LrTDDFT

def Omega(omega_p,n_p,K):
	n = len(omega_p)
	return (omega_p**2*np.identity(n) + 
		2*K*np.sqrt(n_p*omega_p).dot(np.sqrt(n_p*omega_p)))
	

if __name__=='__main__':
	calc = LrTDDFT('lr_dE.dat.gz')
	dump_data(calc,'dump.npz')
	#dump = np:.load('dump.npz')	
