from gpaw.tddft import *

def calc_photoabs_with_kick(direction, kick_strength):
	time_step = 30                  # 1 attoseconds = 0.041341 autime
	total_time = 45 		# femto seconds
	iterations = 45*10**(3)//30

	# Read ground state
	td_calc = TDDFT('../Task1/0bands.gpw')

	# Kick with a delta pulse to z-direction
	td_calc.absorption_kick(kick_strength=kick_strength)

	# Propagate, save the time-dependent dipole moment to 'be_dm.dat',
	# and use 'be_td.gpw' as restart file
	td_calc.propagate(time_step, iterations, 'be_dm_' + direction + '.dat', 
			  'be_td_' + direction + '.gpw')

	# Calculate photoabsorption spectrum and write it to 'be_spectrum_z.dat'
	photoabsorption_spectrum('be_dm_' + direction + '.dat', 'be_spectrum_z.dat')


if __name__=='__main__':
	kick = 10**(-5)
	calc_photoabs_with_kick('x',[kick,0,0])
	calc_photoabs_with_kick('y',[0,kick,0])
