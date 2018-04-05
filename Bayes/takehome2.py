#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import random

#compute m draws from the dirichlet distribition with data given in vector a
def rdirich(m, a):
    th = np.zeros((m, len(a)+1))

    #generates m random gamma draws based on each of the data points
    for j in range(len(a)):
        th[:,j] = random.gamma(a[j], 1, m)
    th = th/np.sum(th, axis=1)[:,np.newaxis] #normalizes
    th[:,-1] = th[:,0] - th[:,1] #calculates gamma
    return th

if __name__ == '__main__':
    #define parameters
    p = 3
    m = int(1e6)

    #define prior and data, calculate posterior
    prior = np.ones(p) * .01
    data = np.asarray((727, 583, 137))
    post = rdirich(m, prior+data)

    #calculate confidence interval size and print with mean
    ci =  np.percentile(post[:,-1], (2.5, 97.5))
    print(np.average(post[:,-1]), ci, np.sum(post[:,-1]>0)/m)

    #plot posterior histograms
    plt.figure(figsize=(10,4))
    plt.subplot(121)
    plt.hist(post[:,0], label = r'$\theta_1$', bins = 50)
    plt.hist(post[:,1], label = r'$\theta_2$', bins = 50)
    plt.xlabel('Vote Fraction')
    plt.ylabel('Counts')
    plt.legend()
    plt.subplot(122)
    plt.hist(post[:,-1], label = r'$\gamma$', bins = 50, color = 'g')
    plt.axvline(x=ci[0], color = 'k', ls = '--')
    plt.axvline(x=ci[1], color = 'k', ls = '--')
    plt.xlabel('Vote Fraction')
    plt.ylabel('Counts')
    plt.legend()
    plt.tight_layout()
    plt.show()
