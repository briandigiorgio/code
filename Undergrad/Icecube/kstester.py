#!/usr/bin/env python

import numpy as np
import argparse
from scipy import special as sp
import time
import matplotlib.pyplot as plt
from scipy import stats as sp

def kstest(d1,d2):
    d1sum = np.sum(d1)
    d2sum = np.sum(d2)
    length = len(d1)
    print len(d1), len(d2)
    if len(d1) != len(d2):
        return "Please input data of the same length"

    cd1 = np.zeros(length)
    cd2 = np.zeros(length)
    c1 = 0
    c2 = 0

    for i in range(int(length)):
        c1 += d1[i]
        c2 += d2[i]
        cd1[i] = c1
        cd2[i] = c2

    #plt.plot(range(len(cd1)),cd1)
    #plt.plot(range(len(cd2)),cd2)
    #plt.show()
    diff = np.absolute(cd1-cd2)
    imax = np.argmax(diff)
    k = diff[imax]
    kn = k*np.sqrt(d1sum)
    p = np.exp((-2*(kn**2)*(d1sum + d2sum)/(d1sum*d2sum)))
    #pstuff = [np.exp(-.125*(((2*i-1)*np.pi)/k)**2) for i in range(100)]
    #p = (np.sqrt(np.pi)/k)*np.sum(pstuff)
    #pstuff = [(-1)**(n-1)*np.exp(-2*(n*k)**2) for n in range(1,100)]
    #p = 2*np.sum(pstuff)
    
    return [p,k,imax]
    
if __name__ == "__main__":
    r1 = np.random.normal(0,1,10000)
    r2 = np.random.normal(2,1,10000)    
    bins = np.arange(-5,5,.01)
    d1 = np.histogram(r1,bins)[0]
    d2 = np.histogram(r2,bins)[0]
    d1 = d1/float(np.sum(d1))
    d2 = d2/float(np.sum(d2))

    #plt.hist(r1)
    #plt.show()
    print len(d1),len(d2)
    print kstest(d1,d2,bins)
    print sp.ks_2samp(r1,r2)
