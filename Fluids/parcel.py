#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import time


#define equation of state
def eos(e, V, g):
    return (g-1) * e / V

#see histogram of data and printout of first 10 elements
def inspect(array):
    plt.hist(array[0], bins = 50, range = (0,L), histtype = 'step', 
            label='x')
    '''
    plt.hist(array[1], bins = 50, range = (0,L), histtype = 'step',
            label='u') 
    plt.hist(array[2], bins = 50, range = (0,L), histtype = 'step',
            label='e')
    plt.hist(array[3], bins = 50, range = (0,L), histtype = 'step',
            label='p')
    plt.hist(array[4], bins = 50, range = (0,L), histtype = 'step',
            label='V')
    plt.hist(array[6], bins = 50, range = (0,L), histtype = 'step',
        label='q')
    '''
    #print(time.time() - start)
    plt.legend()
    np.set_printoptions(suppress = True, threshold = np.nan)
    print("     x   |     u   |    e   |    p   |    V   |    r0  |    q")
    print("-----------------------------------------------------------------")
    print(np.around(array[:, 745:754], decimals = 2).T)
    print(np.nanmax(array[0]))
    plt.show()

#performs a timestep on the input array, outputs the resulting array
#input array of shape 7 by n with indices following format:
#0 = x, 1 = u, 2 = e, 3 = p, 4 = V, 5 = rho0, 6 = q
def timestep(t0, dt):
    #creates new array, carries over initial density
    t1 = t0 
    t1[5] = t0[5]
    t0[2:6,0] = t0[2:6,1]

    #reassigns velocity and position according to equations on handout
    for i in range(n-1):
        t1[1][i] = t0[1][i] - ((2 * dt)/(t0[5][i] + t0[5][i+1])) * \
          (t0[3][i+1] - t0[3][i] + t0[6][i+1] - t0[6][i])/dx

    #dt = np.min(b * dx/(np.abs(t1[1]) + cs))
    print(dt)
    for i in range(n-1):
        t1[0][i] = t0[0][i] + t1[1][i]*dt

        #boundary conditions: rigid wall at 0
        if t1[0][i] <= 0:
            t1[0][i] = 0.1
            t1[1][i] = 0
    
    #recalculate V
    for i in range(n-1):
        t1[4][i] = np.abs((1/t0[5][i]) * (t1[0][i+1] - t1[0][i])/dx)

    #recalculate q for appropriate values
    for i in range(n-1):
        if t1[1][i+1] - t1[1][i] < 0:
            t1[6][i] = .5 * a**2 * (t1[1][i+1] - t1[1][i])**2 * \
                ((1/t1[4][i]) + (1/t0[4][i]))
        else:
            t1[6][i] = 0

    #recalculate e
    for i in range(n-1):
        t1[2][i] = t0[2][i] - (t0[3][i] + t0[6][i]) * \
                (t1[4][i] - t0[4][i]) 

    #recalculate p
    for i in range(n-1):
        t1[3][i] = eos(t1[2][i], t1[4][i], g)

    return t1, dt
if __name__ == "__main__":

    #set parameters for  the simulation
    start = time.time()
    dx = 1
    cs = 1
    b = .01
    dt =  b * dx/cs
    L = 10000
    n = 1000
    a = .5
    g = 7/5

    #set initial conditions: barrier in the middle of the box
    barrier = L//2
    p0 = 1
    r0 = 1
    pl = 1.
    rl = 1.
    pr = .1
    rr = .125

    #create array of particles with characteristics:
    #0 = x, 1 = u, 2 = e, 3 = p, 4 = V, 5 = rho0, 6 = q
    array = np.zeros((7, n))
    
    '''
    #set initial conditions of the particles: all on left of barrier
    j0s = np.linspace(0,barrier, n+1)
    for i in range(n):
        array[0][i] = j0s[i]
        array[3][i] = p0
        array[4][i] = 1/r0
        array[5][i] = r0
    '''

    #j0s = np.concatenate((np.linspace(0,barrier,(3 * n)//4),
    #        np.linspace(barrier+1,L,n//4)))
    j0s = np.linspace(0,L,n)
    for i in range((3 * n)//4):
        array[3][i] = pl
        array[4][i] = 1/rl
        array[5][i] = rl

    for i in range((3*n)//4,n):
        array[3][i] = pr
        array[4][i] = 1/rr
        array[5][i] = rr

    array[0] = j0s
    array[2] = array[3] * array[4]/(g-1)

    #run simulation
    for i in range(10000):
        if not i%10:
            inspect(array)
        array, dt = timestep(array, dt)
        
