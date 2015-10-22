"""
Main script to generate the color-color tracks
"""
#Import all the stuff!!!
import numpy as np
import collections


import matplotlib.pylab as pl
import seaborn as sns; sns.set_style('ticks')
cmap = sns.cubehelix_palette(start=0.0, rot=0.0, light=0.75, dark=0.25, as_cmap=True)
# cmap = sns.color_palette("cubehelix", 7)

from util import synth_mag, synth_mag_pysysp

from specutils.extinction import reddening
from astropy import units as u






def main():
	#Generate data arrays	
	composite = np.genfromtxt('../data/Selsing2015.dat')
	# composite = np.genfromtxt('../data/CompoM.dat')	
	wl, flux, error = composite[:,0], composite[:,1], composite[:,1]
	# dz = np.arange(1.2 , 2.7, 0.05)
	dz = np.arange(0 , 3, 0.01)

	size_redshift = len(dz)
	
	# bands = ['u', 'g', 'r', 'i', 'z', 'Z2', 'Y', 'J', 'H', 'K']
	bands = ['g', 'r', 'J', 'K']


	AV = np.array([0.0, 0.4, 0.8, 1.2, 1.6])
	# AV = np.array([0.0, 0.5, 1.0])
	# AV = np.array([0.0])


	fig, ax = pl.subplots()

	for c, av in enumerate(AV):
		cmap = sns.cubehelix_palette(start=0.5 * c, rot=0.0, light=0.75, dark=0.25, as_cmap=True)
		mags = collections.defaultdict(list)

		#Loop through redshifts
		for ii in dz:
			#Redden at host
			flux_dered =  flux / reddening(wl* u.angstrom, av, r_v=2.72, model='gcc09')
			# error_dered = error / reddening(wl* u.angstrom, av, r_v=2.72, model='gcc09')			
			error_dered= 0
			#Move with redshift
			wl_z = wl * (1+ii)

			#Loop through bands
			for k in bands:
				# mag = synth_mag(band=k, datapath='../data/filter_curves/', wave=wl_z, flux=flux_dered, error=error_dered)
				try:
					mag = synth_mag_pysysp(wl_z, flux_dered, error_dered, bandpath='../data/filter_curves/', band=k)
				except ValueError:
					print('Filter outside composite range for band: '+str(k)+' for redshift' + str(ii))
					mags[k].append(np.nan)
					continue
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
		ax.plot(JK, gJ, label= r'A$_{V}$ = '+str(av), lw=0.3, zorder=1, color=cmap(0.5), alpha= 1)

		sca = ax.scatter(JK, gJ, c=dz, zorder=2, edgecolor='none', alpha = 1, s=5, cmap=cmap)
		print(r'Plotted track for: A$_{V}$ = '+str(av))
	# cmap = sns.cubehelix_palette(start=1, rot=0.0, light=0.75, dark=0.25, as_cmap=True)
	# sca = ax.scatter(JK, gJ, c=dz, zorder=2, edgecolor='none', alpha = 1, cmap=cmap)

	ax.set_xlim((-1.0, 2.0))
	ax.set_ylim((0.0, 4.0))

	cbar = fig.colorbar(sca, ax=ax)
	cbar.set_label('Redshift')

	ax.set_xlabel('J - K')
	ax.set_ylabel('g - J')
	# cbar.set_label('z')
	# set the linewidth of each legend object
	leg = ax.legend(loc=2)
	for legobj in leg.legendHandles:
		legobj.set_linewidth(2.0)

	pl.savefig('../figs/color_track_JKgJ.pdf')
	pl.show()


if __name__ == '__main__':
	main()