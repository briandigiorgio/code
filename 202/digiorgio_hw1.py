#!/usr/bin/env python

###Astro 202 HW1 #3
###Brian DiGiorgio, 2/1/17
###Run this to print out all of the answers

import time
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.optimize import fmin
import warnings 
warnings.filterwarnings("ignore") #only use at the end after debugging

h = 6.62607004081e-27
c = 29979245800
k = 1.3806485279e-16
sigma = (2*np.pi**5*k**4)/(15*c**2*h**3)
T = 10000 #set to whatever


def Bv(nu, T): #Planck equation in terms of frequency
    return ((2*h*(nu**3))/(c**2))/(np.expm1((h*nu)/(k*T)))

def sb(T): #Stefan Boltzmann law
    return sigma*(T**4)

def Bx(x, T): #Planck equation in terms of x = (h*nu)/(kT)
    Cx = (2*(k*T)**3)/((h*c)**2) #Constants out front
    return Cx*x**3/np.expm1(x)

def By(y, T): #Planck equation in terms of y = (hc/(lambda*kT))
    Cy = (2*(k*T)**5)/(h**4*c**3) #Constants out front
    return Cy*y**5/np.expm1(y)

def dBxdT(x, T): #derivative of Bv wrt T translated into x
    #Cx = (2*(k*T)**3)/((h*c)**2)
    #return Cx*(((3*x**2)/np.expm1(x))-((x**3*np.exp(x))/(np.expm1(x)**2)))
    Cx = (2*k**3*T**2)/(h**2)
    return Cx * (x**4 * np.exp(x))/(np.expm1(x)**2)

def part_a(T):
    #perform \int cos\theta I_\nu d\nu d\Omega
    # = \pi \int_0^\inf B_\nu d\nu
    #some corners are cut with integration limits because of errors
    return np.pi * quad(lambda nu: Bv(nu, T), 1e11, 1e16)[0]

def part_b(T):
    #peak values for Planck in terms of dimensionless vars for nu and lambda
    #have to find min values of negatives because of how scipy.optimize works
    xmax = fmin(lambda x: -Bx(x, T), 0, disp = False)[0]
    ymax = fmin(lambda y: -By(y, T), 0, disp = False)[0]
    return xmax, ymax

def part_c(T):
    #find max of derivative in terms of x
    return fmin(lambda x: -dBxdT(x, T), 0, disp = False)[0]

if __name__=="__main__":
    print("Part a:")
    print("Numerical: ", part_a(T))
    print("Stefan-Boltzmann: ", sb(T))
    print("Ratio: ", (part_a(T)/sb(T)), "(definitely close enough)", '\n')

    '''
    x = np.logspace(11,16,100)
    plt.semilogx(x, Bv(x, T))
    plt.show()

    xs = np.logspace(-2,2,100)
    plt.semilogx(xs, By(xs, T))
    plt.show()
    '''

    print("Part b:")
    xmax, ymax = part_b(T)
    print("Peak in x: ", xmax)
    print("Peak in y: ", ymax, '\n')

    print("Part c:")
    print("Peak of derivative: ", part_c(T))
