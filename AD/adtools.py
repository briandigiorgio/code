import time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors, rc, cm
import argparse
import scipy.stats as stats
from scipy.optimize import curve_fit
from astropy import coordinates
import astropy.units as u

def load_fits(file, i = 1):
    return fits.open(file)[i].data.copy()

def loadfits(file, i=1):
    return load_fits(file, i)

def lex(file, i=1):
    hdu = loadfits(file, i=i)
    print('\n'.join(hdu.columns.names))
    return hdu
#reads fits file of data from Kyle, spits out relevant data arrays
#should be of the form 'SPX-GAU-MILESHC-composite_*.fits'
#largely borrowed from Kyle
def read_ad_mass(file, adv = False):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mag = 'elpetro_absmag'
    mage = 'elpetro_abmerr'

    #pull Mi data
    Mi = data[mag][:,5]
    Mie = data[mage][:,5]
    
    #divide by gas rotation speed if user specified
    if adv:
        ad = data['ad2_em']/data['harc_em']
        ade = np.sqrt((data['ad2_se']/data['harc_em'])**2 +
              ((data['ad2_se']*data['harc_se'])/(data['harc_em']**2))**2)

    else: 
        ad = data['ad2_em']
        ade = data['ad2_se']
    
    return Mi, Mie, ad, ade

#reads info on magnitude, ad, color, rotation speed, and identifying info
#taken mostly from Kyle, works on SPX... files
def read_ad_mag(file):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mag = 'elpetro_absmag'
    mage = 'elpetro_abmerr'
    mass = 'elpetro_mass'
    
    #pull Mi data
    Mi = data[mag][:,5]
    Mie = data[mage][:,5]

    #calculate N-r
    Nmr = data[mag][:,1] - data[mag][:,4]
    Nmre = np.sqrt(np.square(data[mage][:,1]) + np.square(data[mage][:,4]))

    #pull rotational velocity data
    gasrc = data['harc_em']
    gasrce = data['harc_se']

    ad = data['ad2_em']
    ade = data['ad2_se']

    plate = data['plate'].astype(str)
    ifu = data['ifudesign'].astype(str)

    return plate, ifu, Mi, Mie, Nmr, Nmre, gasrc, gasrce, ad, ade

def find_repeats(array):
    return np.setdiff1d(array, np.unique(array))

def powerlaw(x, a, b, c):
    return a * np.power(x, b) + c

def line(x, m, b):
    return m * x + b

def exponential(x, a, b):
    return a * b ** x

def weighted_avg_std(values, weights):
    average = np.average(values, weights=weights)
    sumweights = np.sum(weights)
    variance = np.sum((weights * (values - average)**2))/sumweights
    return (average, np.sqrt(variance)/np.sqrt(len(values)))

#makes an array of length numlines going through the specified color map
#plt.plot(..., c = colors[i]) when plotting multiple lines
def make_cmap(numlines, cmap):
    cnorm = colors.Normalize(vmin = 0, vmax = numlines)
    scalarmap = cm.ScalarMappable(norm = cnorm, cmap = cmap)
    return scalarmap.to_rgba(range(numlines))

#returns an array of indices of cat2 that correspond to the same item in cat1
def match_cats(ra1, dec1, ra2, dec2):
    cat1 = coordinates.SkyCoord(ra = ra1 * u.degree, dec = dec1 * u.degree)
    cat2 = coordinates.SkyCoord(ra = ra2 * u.degree, dec = dec2 * u.degree)
    #cat1tocat2, d2d, d3d = cat1.match_to_catalog_sky(cat2)
    cat1tocat2, d2d, d3d = coordinates.match_coordinates_sky(cat1, cat2)
    '''
    print(len(cat1tocat2))
    bad = []
    for i in range(len(d2d)):
        if d2d[i].is_within_bounds('10m',None):
            bad += [i]
    cat1tocat2 = np.delete(cat1tocat2, bad)
    print(len(cat1tocat2), len(bad))
    '''
    return cat1tocat2, d2d

#reads fits file of data from Kyle, spits out relevant data arrays
#should be of the form 'SPX-GAU-MILESHC-composite_*.fits'
#largely borrowed from Kyle
def read_mass(file, adv = False):

    #loads file
    data = load_fits(file)

    #names of columns to be pulled
    mass = 'elpetro_mass'
    #mass = 'sersic_mass'

    #pull Mi data
    m = data[mass]
    
    #divide by gas rotation speed if user specified
    if adv:
        ad = data['ad2_em']/data['harc_em']
        ade = np.sqrt((data['ad2_se']/data['harc_em'])**2 +
              ((data['ad2_se']*data['harc_se'])/(data['harc_em']**2))**2)

    else: 
        ad = data['ad2_em']
        ade = data['ad2_se']
    
    return m, ad, ade


