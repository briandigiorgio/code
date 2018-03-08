#!/usr/bin/env python

###THIS IS BASICALLY A TESTBED AT THIS POINT
###TO ACTUALLY GET THE ENERGY FOR A PIXEL, USE POWERLAW.PY

import numpy as np
import matplotlib.pyplot as plt
import argparse
import scipy.optimize as sp
from scipy import stats

desc = "takes one pixel from a histogram and plots its NChannel distribution"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str, help = ".npy file to read")
parser.add_argument("-p", "--pixel", type = int, help = "pixel number")
parser.add_argument("-l", "--lower", type = int, default = 0, 
                    help = "lower cut of bins to find fitline of")
parser.add_argument("-u", "--upper", type = int, default = -1,
                    help = "upper cut of bins to find fitline of")
args = parser.parse_args()

#curve for fitting
def curve(x, a, b):
    return a * x ** b

if __name__ == "__main__":

    #load file, define variables for fit
    matrix = np.load(args.file)
    data = matrix[args.pixel]

    for i in range(len(data)):
        if data[i] == 0:
            data[i] = .00001

    sum = float(np.sum(data))
    normed = data/sum
    error = np.sqrt(data)/sum

    #give 0 points a large error to make them not matter
    for i in range(len(error)):
        if error[i] == 0.:
            error[i] = np.inf

    #generate modified exponential y bins
    print "Generating with modified exponential y axis..."
    ny = 100
    yrange = np.arange(ny-2)
    ybins = [6.,7.]
    offset = 0.
    step = 1.

    #generate 100 bins ranging from 0 to 1000 that grow exponentially
    for i in yrange:
        newbin = np.trunc(np.exp(np.log(1000) * i/(ny-2))) + offset
        ybins.append(newbin)
        diff = ybins[-1] - ybins[-2]

        #if the size of the bin has grown, make that the min bin size
        if diff > step:
            step = diff

        #increase the bin size until it is the minimum bin size
        while diff < step:
            ybins[-1] += 1
            offset += 1
            diff += 1

    #calculate binwidth for exponential y bins 
    binwidth = np.zeros(len(ybins))
    for i in range(len(binwidth)):
        if i+1 > len(ybins)-1:
            binwidth[i] = np.inf
        else:
            binwidth[i] = ybins[i+1] - ybins[i]

    #reset y bins and bin width if using true y axis
    if matrix.shape[1] == 1000:
        ybins = np.linspace(0,1000,1000)
        binwidth = np.zeros(1000) + 1

    #normalize data with binwidths, plot it with error bars
    normed = normed / binwidth
    if args.upper ==-1:
        args.upper = 1000
    lower,upper=np.digitize([args.lower,args.upper],ybins)

    print "Plotting..."
    plt.loglog(ybins, normed, "b.")
    plt.errorbar(ybins, normed, yerr = error, fmt = '.')

    #calculate best fit power law, plot it
    x = np.linspace(0,1000,1000)

    popt, pcov = sp.curve_fit(curve,ybins[lower:upper],normed[lower:upper],
                              sigma=error[lower:upper])
    plt.loglog(x, curve(x, popt[0], popt[1]))
    #plt.errorbar(ybins, normed, yerr = error/binwidth, fmt = '.')
    plt.xlim(args.lower,args.upper)
    plt.ylim(ymax = 1)
    print popt, pcov
    print "Estimated curve: %fx^%f" % (popt[0],popt[1])
    print "Standard Deviations: " + str(np.sqrt(np.diag(pcov)))


    for i in range(len(normed)):
        if normed[i] == 0:
            normed[i] = 1
            print "replaced"
    '''
    logbins = np.log10(ybins)
    lognormed = np.log10(normed)
    plt.plot(logbins, lognormed, "b.")


    if not args.upper:
        args.upper = len(logbins) - 1
    lower, upper = np.digitize([np.log10(args.lower), np.log10(args.upper)], logbins)

    logbins = logbins[lower:upper]
    lognormed = lognormed[lower:upper]
    #scalederror = (error/normed) * lognormed
    #plt.errorbar(logbins, lognormed, yerr = scalederror, fmt = '.')
    popt, pcov = sp.curve_fit(curve,logbins,lognormed)#,sigma = error)
    print popt
    x = np.linspace(0,3,1000)
    plt.plot(x, curve(x, popt[0], popt[1]))
    error = np.zeros((2,len(data)))
    error[0] = (data + np.sqrt(data))/data
    error[1] = (data - np.sqrt(data))/data
    print error

    scalederror = np.zeros(error.shape)
    scalederror = [normed*(1 - error[0]), normed*(1 + error[1])]
    logerror = np.log10(scalederror)
    #plt.errorbar(logbins, lognormed, yerr = logerror, fmt = ",")
    '''

    plt.show()
