"""
Functions for the calculation of color-color tracks
"""

__all__ = ["synth_mag"]

def synth_mag(band=None, datapath=None, wave=None, flux=None):
    """
	Function to calculate synthetic AB magnitudes. Following Bessell & Murphy (2012) and Casagrande & VandenBerg (2014). 


    Args:
        band (string): string value for the band. Must match .dat files for the transmission curves
        datapath (string): Path to transmission curves
        wave (np.array): Input wavelength array
        flux (np.array) = Input flux array

    Returns:
        float: Calculated apparent magnitude

    """

    import glob
    import numpy as np
    from astropy.cosmology import Planck13 as cosmo
    from gen_methods import medfilt

    filter = glob.glob(datapath+band+'.dat')[0]
    filter = np.genfromtxt(filter)
    wl_filt = np.array(filter[:,0])
    if np.mean(wl_filt) < 10:
        wl_filt *= 10000.0
    filt = filter[:,1]
    filt_new =  common_wavelength(wl_filt, wave, filt, fill_value=0.0)
    prod = filt_new * flux
    numerator = np.sum(prod * wave)
    denom = 3e18 * np.sum(filt_new/wave)
    f_nu = numerator / denom
    i_band_mag = -2.5 * np.log10(f_nu) - 48.6
    return i_band_mag



def common_wavelength(wlarr_old, wlarr_new, fluxarr_old, fill_value = 0.):
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