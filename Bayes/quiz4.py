#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammaln

#define global variables
n = 672
thh = 191.8
x = np.linspace(160,240,81)

#log likelihood
def logl(th, n, thh):
    return (-n/2)*np.log(th) - n*thh/(2*th)

#make plot of log likelihood
def num1():
    plt.plot(x, logl(x, n, thh))
    plt.xlabel('Variance')
    plt.ylabel('Log Likelihood')
    plt.show()

#scaled inverse chi square density, lifted from David's online R code
def invchisq(th, nu, th0):
    return np.exp((nu/2)*np.log(nu/2) - gammaln(nu/2) + (nu/2)*np.log(th0) - \
            (1+nu/2)*np.log(th) - nu*th0/(2*th))

#plot prior, likelihood, and posterior
def num3():
    plt.plot(x, 2.5/x, label = 'Prior')
    plt.plot(x, invchisq(x, n-2, (n/(n-2))*thh), label = 'Likelihood')
    plt.plot(x, invchisq(x, n, thh), label = 'Posterior')
    plt.xlabel('Variance')
    plt.ylabel('Likelihood')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    num1()
    num3()
