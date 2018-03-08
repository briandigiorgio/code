#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import random

def rdirich(m, a):
    if len(a) != 3:
        print('nope, only handles p=3')
        return
    
    th = np.zeros((m, len(a)))

    th[:,0] = random.gamma(a[0], 1, m)
    th[:,1] = random.gamma(a[1], 1, m)
    th[:,2] = random.gamma(a[2], 1, m)

    th[:,0] = th[:,0]/np.sum(th[:,0])
    th[:,1] = th[:,1]/np.sum(th[:,1])
    th[:,2] = th[:,2]/np.sum(th[:,2])
    return th

if __name__ == '__main__':
    p = 3
    prior = np.ones(p) * .01
    data = np.asarray((727, 583, 137))
    print(rdirich(10000, prior+data)) 
