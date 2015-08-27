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
# cmap = sns.cubehelix_palette(n_colors=6, start=0.0, rot=1.5, gamma=1.0, hue=1.0, light=0.85, dark=0.15, reverse=True, as_cmap=False)
cmap = sns.color_palette("cubehelix", 3)

from util import synth_mag



def main():
	#Generate data arrays	
	composite = np.genfromtxt('../data/Selsing2015.dat')
	wl, flux, error = composite[:,0], composite[:,1], composite[:,2]
	dz = np.arange(0, 5, 0.01)
	bands = ['u', 'g', 'r', 'i', 'z', 'Z2', 'Y', 'J', 'H', 'K']
	mags = collections.defaultdict(list)

	#Loop through redshifts
	for i in dz:
		wl_z = wl * (1+i)
		flux_z = flux / (1+i)
		error_z = error / (1+i)

		#Loop through bands
		for k in bands:
			mag = synth_mag(band=k, datapath='../data/filter_curves/', wave=wl_z, flux=flux_z)
			mags[k].append(mag)


	#Calculate colors
	JK = np.array(mags['J']) - np.array(mags['K'])
	gJ =  np.array(mags['g']) - np.array(mags['J'])


	#Printing
	print('J - K:', JK)
	print('J:', np.array(mags['J']))
	print('g - J:', gJ)

	#Make the plot
	fig, ax = pl.subplots()
	sca = ax.scatter(JK, gJ, c=dz)
	ax.set_xlim((-1, 2))
	ax.set_ylim((-0.5, 2.0))
	cbar = fig.colorbar(sca, ax=ax)
	ax.set_xlabel('J - K')
	ax.set_ylabel('g - J')
	cbar.set_label('z')
	pl.savefig('../figs/color_track.pdf')
	pl.show()


if __name__ == '__main__':
	main()