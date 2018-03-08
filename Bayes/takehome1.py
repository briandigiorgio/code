#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.special import gamma
from scipy.integrate import quad

def num2iia(e,pi):
    p = np.linspace(0,.1,100)
    ppv = p*e/(p*e + (1-p)*(1-pi))
    npv = (1-p)*pi/(p*(1-e) + (1-p)*pi)
    plt.plot(p, ppv, label = 'PPV')
    plt.plot(p, npv, label = 'NPV')
    plt.legend()
    plt.grid(True)
    plt.show()

def ppv(p, e, pi):
    return p*e/(p*e + (1-p)*(1-pi))

def npv(p, e, pi):
    return (1-p)*pi/(p*(1-e) + (1-p)*pi)

def num2iii():
    print("a. ", ppv(.01, 1, .98))
    print("b. ", ppv(.01, .95, .997))
    print("c. ", npv(.01, 1, .98), npv(.01, .95, .997))

def lexp(l, n, yb):
    return (1/(l**n))*np.exp((-n*yb)/l)
    #return (l**n)*np.exp(l*n*yb)

def llexp(l, n, yb):
    return ((-n*yb)/l) - n*np.log(l)

def cdfi(i, n, yb):
    return -yb*np.log(1-(i-.5)/n)

def numBib():
    y = (495,541,1461,1555,1603,2201,2750,3468,3516,4319,6622,7728,13159,21194)
    yb = np.average(y)
    n = len(y)
    print(yb)
    l=np.linspace(2000,15000,100)
    plt.subplot(211)
    plt.plot(l, lexp(l, n, yb))
    plt.axvline(yb, ls = '--', c='red', label = 'Mean')
    plt.title('Likelihood')
    plt.legend()

    plt.subplot(212)
    plt.plot(l, llexp(l, n, yb))
    plt.axvline(yb, ls = '--', c='red')
    plt.title('Log Likelihood')

    plt.tight_layout()
    plt.show()

    x = np.linspace(0,12000,100)
    plt.plot(cdfi(np.asarray(list(range(n))), n, yb), y, 'o')
    plt.plot(x,x,'k--')
    plt.title('Exponential Probability Plot')
    ax = plt.gca()
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()

def igamma(l, a, b):
    return (b**a)/gamma(a) * l**(-a-1)*np.exp(-b/l)

def numBvc():
    x = np.linspace(1000,12000,100)
    plt.plot(x, igamma(x, 8.25, 32625), label = 'Prior')
    plt.plot(x, igamma(x, 22.25, 103237), label = 'Posterior')
    plt.plot(x, igamma(x, 13, 70612), label = 'Likelihood')
    plt.legend()
    plt.xlabel('Lambda')
    plt.ylabel('Probability')
    ax = plt.gca()
    ax.axes.get_yaxis().set_ticks([])
    plt.show()


if __name__ == '__main__':
    #num2iia(.95, .98)
    #num2iii()
    #numBib()
    numBvc()
