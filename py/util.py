"""
Functions for the calculation of color-color tracks
"""

__all__ = ["synth_mag"]

def synth_mag(band=None, datapath=None, wave=None, flux=None):
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
    # prod = medfilt(filt_new * flux, 15)
    prod = filt_new * flux
    numerator = np.sum(prod * wave)
    denom = np.sum(filt_new * (3e18/wave))
    f_nu = numerator / denom
    i_band_mag = -2.5 * np.log10(f_nu) - 48.6
    return i_band_mag


def common_wavelength(wlarr_old, wlarr_new, fluxarr_old, fill_value = 0.):
    from scipy import interpolate
    f = interpolate.interp1d(wlarr_old, fluxarr_old, kind='linear', bounds_error = False, fill_value=fill_value)
    fluxarr_new = f(wlarr_new)
    return fluxarr_new