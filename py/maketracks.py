"""
Main script to generate the color-color tracks
"""
#Import all the stuff!!!
import numpy as np
import collections

from matplotlib import rc_file
rc_file('/Users/jselsing/Pythonlibs/plotting/matplotlibstyle.rc')
import matplotlib.pylab as pl
import seaborn as sns; sns.set_style('ticks')
cmap = sns.cubehelix_palette(n_colors=7, start=0.0, rot=1.5, gamma=1.0, hue=1.0, light=0.85, dark=0.15, reverse=True, as_cmap=False)
# cmap = sns.color_palette("cubehelix", 7)

from util import synth_mag

from specutils.extinction import reddening
from astropy import units as u




def main():
	#Generate data arrays	
	composite = np.genfromtxt('../data/Selsing2015.dat')
	wl, flux, error = composite[:,0], composite[:,1], composite[:,2]
	dz = np.arange(0, 5, 0.01)
	# bands = ['u', 'g', 'r', 'i', 'z', 'Z2', 'Y', 'J', 'H', 'K']
	bands = ['g', 'r', 'J', 'K']


	AV = np.array([0.0, 0.4, 0.8, 1.2, 1.6, 2.0])
	# AV = np.array([0.0])
	# EBV = AV
	# colors = 
	fig, ax = pl.subplots()

	for c, av in enumerate(AV):
		mags = collections.defaultdict(list)

		#Loop through redshifts
		for i in dz:
			#Redden at host
			flux_dered =  flux / reddening(wl* u.angstrom, av, r_v=2.72, model='gcc09')
			error_dered = error / reddening(wl* u.angstrom, av, r_v=2.72, model='gcc09')			

			wl_z = wl * (1+i)

			#Loop through bands
			for k in bands:
				mag = synth_mag(band=k, datapath='../data/filter_curves/', wave=wl_z, flux=flux_dered)
				mags[k].append(mag)


		#Calculate colors
		JK = np.array(mags['J']) - np.array(mags['K'])
		gJ =  np.array(mags['g']) - np.array(mags['J'])
		gr = np.array(mags['g']) - np.array(mags['r'])



		#Printing
		# print('J - K:', JK)
		# print('J:', np.array(mags['J']))
		# print('g - J:', gJ)

		#Make the plot
		ax.plot(JK, gJ, color=cmap[c], label= r'A$_{V}$ = '+str(av), lw=2.0, zorder=1)	
		sca = ax.scatter(JK, gJ, c=dz, zorder=2, alpha=0.5)
		print(r'Plotted track for: A$_{V}$ = '+str(av))

	ax.set_xlim((-1.0, 2.0))
	ax.set_ylim((0.0, 4.0))
	cbar = fig.colorbar(sca, ax=ax)
	ax.set_xlabel('J - K')
	ax.set_ylabel('g - J')
	cbar.set_label('z')
	pl.legend(loc=2)
	pl.savefig('../figs/color_track_JKgJ.pdf')
	pl.show()


if __name__ == '__main__':
	main()