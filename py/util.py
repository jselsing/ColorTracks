"""
Functions for the calculation of color-color tracks
"""

__all__ = ["synth_mag", "synth_mag_pysysp"]

import numpy as np

def synth_mag(band=None, datapath=None, wave=None, flux=None, error=None):
    """
	Function to calculate synthetic AB magnitudes. Following Bessell & Murphy (2012) and Casagrande & VandenBerg (2014). 


    Args:
        band (string): string value for the band. Must match .dat files for the transmission curves
        datapath (string): Path to transmission curves
        wave (np.array): Input wavelength array
        flux (np.array) = Input flux array
        error (np.array) = Input error array

    Returns:
        float: Calculated apparent magnitude

    """

    import glob
    import numpy as np
    from astropy.cosmology import Planck13 as cosmo
    from gen_methods import medfilt
    from supersmoother import SuperSmoother

    #Read file
    filter = glob.glob(datapath+band+'.dat')[0]
    filter = np.genfromtxt(filter)
    filter_wave = np.array(filter[:,0])
    filter_throughput = filter[:,1]

    #Convert micron to Angstrom
    if np.mean(filter_wave) < 10:
        filter_wave *= 10000.0

    #Interpolate filter to same sampling as target spectrum
    filt_new =  common_wavelength(filter_wave, filter_throughput, wave, fill_value=0.0)

    #Smooth flux to reject noisy pixels:
    sflux = medfilt(flux, 25)

    #Calculate AB magnitude
    f_nu = (np.sum(filt_new * sflux * wave)) / (3e18 * np.sum(filt_new/wave))
    i_band_mag = -2.5 * np.log10((f_nu)) - 48.6
    return i_band_mag



def common_wavelength(wlarr_old, fluxarr_old, wlarr_new, fill_value = 0.):
    """
	Function to interpolate flux arrays onto new grid.


    Args:
        wlarr_old (np.array): Input wavelength array to be interpolated.
        wlarr_new (np.array): Target grid
        fluxarr_old (np.array): Flux to be interpolated into target grid
        fill_value (float) = Fill value for values outside interpolation range

    Returns:
        np.array: Resampled flux array

    """
    from scipy import interpolate
    f = interpolate.interp1d(wlarr_old, fluxarr_old, kind='linear', bounds_error = False, fill_value=fill_value)
    fluxarr_new = f(wlarr_new)
    return fluxarr_new


import sys
sys.path.append('/Users/jselsing/github/pysysp/pysysp/')
import pysysp


def synth_mag_pysysp(wl, flux, error, bandpath, band):
    fil = pysysp.StarSpectrum()
    fil.setwavelength(wl)
    fil.setflux(flux)
    i = pysysp.BandPass(bandpath+band+'.dat')
    if np.mean(i.wavelength) < 10:
        i.wavelength *= 10000
        i.response /= 10000
        i.update_bandpass()
    # print(i.wavelength)
    return fil.apmag(band=i, mag='AB')


    





