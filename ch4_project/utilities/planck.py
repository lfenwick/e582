"""
   calculate the planck function and the brightness temperature.  Note that
   wavelength and radiance are in Modis units (microns, W/m^2/micron/sr) not
   mks.  See http://en.wikipedia.org/wiki/Moderate-Resolution_Imaging_Spectroradiometer
   for modis channel wavelength ranges.
"""


import numpy as np

c=2.99792458e+08  #m/s -- speed of light in vacumn
h=6.62606876e-34  #J s  -- Planck's constant
kb=1.3806503e-23  # J/K  -- Boltzman's constant
c1=2.*h*c**2.
c2=h*c/kb

def planckwavelen(wavel,Temp):
    """input wavelength in microns and Temp in K, output
    bbr in W/m^2/micron/sr
    """
    wavel=wavel*1.e-6  #convert to meters
    Blambda=1.e-6*c1/(wavel**5.*(np.exp(c2/(wavel*Temp)) -1))
    return Blambda

def planckInvert(wavel,Blambda):
    """input wavelength in microns and Blambda in W/m^2/micron/sr, output
    output brightness temperature in K
    """
    Blambda=Blambda*1.e6  #convert to W/m^2/m/sr
    wavel=wavel*1.e-6  #convert wavelength to m
    Tbright=c2/(wavel*np.log(c1/(wavel**5.*Blambda) - 1.))
    return Tbright
     
