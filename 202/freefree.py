#!/usr/bin/env python

import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

#set physical constants (all cgs)
e = 4.80320425e-10 #fundamental charge
h = 6.62606885e-27 #planck
c = 2.99792458e10 #speed of light
me = 9.10938356e-28 #electron mass
kb = 1.38064852e-16 #boltzmann
ry = 2.1798723e-11 #rydberg energy
Z = 1 #nuclear charge, assumed to be hydrogen
n = 1 #density, assumed n_i and n_{e^-} are the same

#return a coordinate grid with ordered pairs of each x,y combination
#equivalent to [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]
#barely even faster, so there's probably a more elegant way to do this
def makegrid(x, y):
    points = np.asarray(np.meshgrid(x, y, indexing = 'ij'))
    return np.swapaxes(np.swapaxes(points, 0, 1), 1, 2)

#calculate the gaunt free free factor from input data table at specified point
def gff(table, u, g2):
    #make coordinate space for table and change to be right shape
    us = np.logspace(-16,13,146)
    g2s = np.logspace(-6,10,81)
    points = makegrid(us, g2s).reshape(len(us)*len(g2s), 2)

    #interpolate between table points in most accurate way to get value
    return griddata(points, table, (u,g2), method = 'cubic')

#return j for given range of nus at some T
def emission(nus, T):
    #calculate values relevant to gaunt
    u = (h*nus)/(kb*T)
    inspect(u)
    g2 = ((Z**2)*ry)/(kb*T) 
    print(g2)

    #read in data table (van Hoof+ 2014) and get gff values 
    table = np.loadtxt('gauntff.dat')[:146,:].flatten()
    gffs = gff(table, u, g2)
    inspect(gffs)
    #implement equation 5.14b from Rybicki & Lightman to get j
    C = ((32*e**4*h)/(12*np.pi*me**2*c**3)) * np.sqrt((np.pi*ry)/(3*kb))#5.44
    return C * n**2 * Z**2 * T**(-1/2) * np.exp((-h*nus)/(kb*T)) * gffs

#calculate thermal free free absorption from Rybicki/Lightman equation 5.18b
def absorption(nu, T):
    u = (h*nus)/(kb*T)
    g2 = ((Z**2)*ry)/(kb*T) 
    
    table = np.loadtxt('gauntff.dat')[:146,:].flatten()
    gffs = gff(table, u, g2)

    C = (4*e**6)/(3*me*h*c) * np.sqrt((2*np.pi)/(3*kb*me))
    return C * T**(-1/2) * Z**2 * n**2 * nus**-3 * (1-np.exp((-h*nus)/(kb*T))) * gffs

def inspect(array):
    print(array[0], array[len(array)//2], array[-1])

if __name__ == '__main__':
    #set temp, set up lambdas, calculate nus
    T = 10000
    lambdas = np.logspace(-1,1,100) * 1e-4
    inspect(lambdas)
    nus = c/lambdas
    inspect(nus)

    #calculate and plot j against lambda
    js = emission(nus,T)
    nugnu = 4*np.pi * nus * js
    plt.loglog(lambdas, nugnu, 'k-')
    #plt.loglog(lambdas, absorption(nus, T))
    #plt.loglog(lambdas, nugnu/absorption(nus, T))
    plt.ylabel(r'$\nu\, \gamma_\nu$')
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.xticks([1e-5, 1e-4, 1e-3], [0.1, 1, 10])
    plt.annotate('T = %d K' % T, xy = (.8, .1), xycoords = 'axes fraction')
    plt.annotate(r'$N_i = N_{e^-} = $%d' % n, xy = (.8, .05), xycoords = 'axes fraction')
    plt.show()
