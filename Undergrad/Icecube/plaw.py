#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
import time
import matplotlib.pyplot as plt
from scipy import optimize as sp

desc = "Takes a .npy 2d histogram output by dataprocessor and calculates" \
     +  "the best fit power law for each pixel, plots them"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str,
                    help = ".npy file to be read")
parser.add_argument("--min", type = float, help = "minimum plot value")
parser.add_argument("--max", type = float, help = "maximum plot value")
parser.add_argument("--moll", default = False, action = "store_true", 
                    help = "plot in mollview instead of orthview")
parser.add_argument("-b", "--bandsize", type = int, default = 0, 
                    help = "the size of the declination bands in pixels, "
                    + "odd numbers make more sense")
parser.add_argument("-p", "--pixel", type = int, default = 0,
                    help = "display histogram for one pixel with fits")
parser.add_argument("-l", "--lower", type = int, default = 0, 
                    help = "lower cut of bins to find fitline of")
parser.add_argument("-u", "--upper", type = int, default = -1,
                    help = "upper cut of bins to find fitline of")
parser.add_argument("--save", type = str, default = "", 
                    help = "saves plot to filepath given instead of showing it")
parser.add_argument("--mask", type = int, default = 0,
                    help = "number of degrees to mask around the pole")
parser.add_argument("--stdev", default = False, action = "store_true",
                    help = "shows plot of stdev instead")
parser.add_argument("--sigma", default = False, action = "store_true",
                    help = "shows plot of how many stdevs away from 0 each"
                    + "pixel is") 
args = parser.parse_args()

#fit line function
def curve(x, a, b):
    return a * x ** b

#function to find a fit line for a given set of data, returns slope
def fit(data, ybins, curve,  bandsize = 0):

    #divides by bandsize if given to renormalize data
    if bandsize:
        data = data/bandsize

    #find total counts in data to normalize data and error
    sum = float(np.sum(data))
    normed = data/sum
    error = np.sqrt(data)/sum

    #give 0 points a large error so they don't matter
    for i in error:
        if not i:
            i = np.inf

    #set upper and lower bounds for cut of fit line, cut the data
    if not args.upper:
        args.upper = len(logbins) - 1

    #find the upper and lower limits in terms of bins, slice data
    lower,upper=np.digitize([args.lower,args.upper],ybins)
    bins = ybins[lower:upper]
    normed = normed[lower:upper]
    error = error[lower:upper]

    #return best fit line slope
    return sp.curve_fit(curve, bins,normed,sigma = error)

if __name__ == "__main__":
    start = time.time()

    print "Importing..."
    matrix = np.load(args.file)

    #define parameters
    NPIX = matrix.shape[0]
    BINS = matrix.shape[1]
    NSIDE = hp.npix2nside(NPIX)

    #generate y bins with modified exponential bins 
    print "Generating bins..."
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

    print "Finding bands..."
    #construct a list of zenith bands by scanning for unique theta values
    thetas = hp.pix2ang(NSIDE, range(NPIX))[0]
    thetatest = thetas[0]
    zbands = np.unique(thetas)

    #find the appropriate bin number for each pixel
    binnum = np.zeros(NPIX)
    band = 0
    for pixel,theta in enumerate(thetas):
        binnum[pixel] = np.digitize([theta], zbands) - 1

    #find width of each bin for normalization purposes
    binwidth = np.zeros(len(ybins))
    for i in range(len(binwidth)):
        if i+1 > len(ybins)-1:
            binwidth[i] = np.inf
        else:
            binwidth[i] = ybins[i+1] - ybins[i]
   
    #make arrays of sums for each bin and the power law exponent in each band
    print "Extracting..."
    zenavg = np.zeros((len(zbands), BINS))
    pavgs = np.zeros(len(zbands))
    varavgs = np.zeros(len(zbands))

    #add the pixels in each bin together for sum, calculate best fit for each
    for i, zenbin in enumerate(zbands):
        zcut = (thetas == zenbin)
        zenavg[i] += np.sum(matrix[zcut], axis=0)/float(np.sum(zcut))

    #hide masked pixels
    dp = np.zeros(NPIX)
    stdev = np.zeros(NPIX)

    #if true bins, reassign bin width and bins
    if BINS == 1000:
        binwidth = np.zeros(matrix.shape[1]) + 1 
        ybins = np.linspace(0,1000,1000)

    #find difference in best fit for each pixel
    print "Calculating..."
    for pixel,zenbin in enumerate(binnum):
        average = np.zeros(len(zenavg[zenbin]))
        average += zenavg[zenbin]

        #only run on pixels that are less than 65 deg from pole
        if thetas[pixel] >= 2:

            #widens the declination bands if value hasn't been calculated yet
            if args.bandsize and not pavgs[zenbin]:
                average = np.zeros(len(zenavg[zenbin]))
                average += zenavg[zenbin]
                size = 1
                delta = 1

                #while the bands are not as big as you want them, widen them
                #adds a row of pixels above then below then above then below...
                while size < args.bandsize:

                    #if band is pushed over pole, bring it back
                    if zenbin + delta >= 255 or zenbin + delta < 0:
                        delta = -delta - 1

                    #actually adds the row avarages together
                    average += zenavg[zenbin + delta]

                    #adjusts the offset 
                    size += 1
                    delta = -delta
                    if delta > 0:
                        delta += 1 

                #calculate and store slope of fit line
                opt,cov = fit(average, ybins, curve, args.bandsize)
                pavgs[zenbin] = opt[1]
                varavgs[zenbin] = cov[1][1]

            #if value hasn't been calculated already, run through same process
            elif not pavgs[zenbin]:
                average = np.zeros(len(zenavg[zenbin]))
                average += zenavg[zenbin]
                pavgs[zenbin] = fit(average, ybins, curve, 0)

            #finds difference in power law exponent for each pixel
            popt,pcov = fit(matrix[pixel], ybins, curve)
            pixfit = popt[1]
            var = pcov[1][1]
            dp[pixel] = pixfit - pavgs[zenbin]
            stdev[pixel] = np.sqrt(var + varavgs[zenbin])


            if not pixel % 1000 and pixel:
                print "Pixels processed: %d,000" % (pixel/1000)
                print dp[pixel], stdev[pixel], dp[pixel]/stdev[pixel]
    print "Time taken: %f seconds" % (time.time() - start)
    
    #hide pixels in the northern hemisphere and plot 
    print "Plotting..."
    
    if args.stdev:
        param = stdev
    if args.sigma:
        param = dp/stdev
    else:
        param = dp

    #hide pixels that are more than 65 deg from pole or are too close to pole
    for i in range(NPIX):
        if thetas[i] < 2:
            dp[i] = hp.UNSEEN

        #mask weirdness at pole
        if args.mask and (np.pi - thetas[i]) * 180/np.pi < args.mask:
            dp[i] = hp.UNSEEN

    #show in IC standard rotated half sky mollview
    if args.moll:
        hp.mollview(param, min = args.min, max = args.max, cbar = False, rot = 180)
        plt.ylim(-1, .005)

    #show an individual pixel's histogram compared to average if specified
    #error bars not currently working
    elif args.pixel:
        
        #get data, fix errors
        data = matrix[args.pixel]

        for i in range(len(data)):
            if data[i] <= 0:
                data[i] = .0001

        #normalize, get error bars (not working yet)
        sum = float(np.sum(data))
        normed = data/sum
        normed = normed/binwidth
        error = np.sqrt(data)/sum


        #recaluculate best fit for band, plot it
        bdata = zenavg[binnum[pixel]]
        bsum = np.sum(bdata)
        bnormed = bdata/bsum
        bnormed = bnormed / binwidth
        berror = np.sqrt(bdata)/bsum

        for i in range(len(normed)):
            if normed[i] == 0:
                normed[i] = 1

            if bnormed[i] == 0:
                bnormed[i] = 1

        #plot everything
        plt.loglog(ybins, normed, "b.")
        plt.plot(ybins, bnormed, "r.")
        plt.errorbar(ybins, normed, yerr = error, fmt = ",")

        #set upper and lower limits and slice
        if not args.upper:
            args.upper = len(logbins) - 1

        lower,upper=np.digitize([np.log10(args.lower),np.log10(args.upper)],logbins)

        logbins = logbins[lower:upper]
        lognormed = lognormed[lower:upper]
        blognormed = blognormed[lower:upper]

        #find fit lines of slices
        popt, pcov = sp.curve_fit(curve,logbins,lognormed)#,sigma = error)
        popt2, pcov2 = sp.curve_fit(curve,logbins,blognormed)
        print popt, popt2

        #plot fitlines
        x = np.linspace(0,3,1000)
        plt.plot(x, curve(x, popt[0], popt[1]), lw = 1.5, color = "blue")
        plt.plot(x, curve(x, popt2[0], popt2[1]), lw = 1.5, color = "red")
        plt.xlim(.75,3)
        plt.ylim((-10, 0))

    #show in normal IceCube orthview of southern hemisphere if not specified
    else:
        hp.orthview(param, rot = [0,-90,180], half_sky = True, min = args.min,
                    max = args.max, cbar = False)

    #gridlines
    hp.graticule()
    '''
    #a whole long thing to generate a caption for the plot based off args
    print
    title = "Power Spectrum fit plot for " + args.file
    if args.smoothing:
        title += " with %d degrees of smoothing" % args.smoothing
    if args.bandsize:
        if args.smoothing:
            title += " and"
        title += " with %d pixel declination bands" % args.bandsize
    if args.merge:
        if args.smoothing or args.bandsize:
            title += " and"
        title += " with NChannel bins of width %d" % args.merge
    if args.lower:
        if args.smoothing or args.bandsize or args.merge:
            title += " and"
        title += " a lower NChannel limit of %d" % args.lower
    if args.upper != -1:
        if args.lower:
            title += " and an upper limit of %d" % args.upper
        elif args.smoothing or args.bandsize or args.merge:
            title += " and an upper Nchannel limit of %d" % args.upper
    if args.pixel:
        " for pixel number %d" % args.pixel
    if args.normed:
        title += ", normalized"
    
    #make each line of the caption the correct length
    for i in range(len(title)/80):
        title = title[:79 * (i+1)] + "- \n" + title[79 * (i+1):]

    #add title and caption
    print title
    #plt.figtext(.1, .02, title)
    '''
    if not args.pixel:
        plt.title("Energy Spectrum Deviation")
     
        #make color bar correct size
        fig = plt.gcf()
        ax = plt.gca()
        image = ax.get_images()[0]
        cbar = fig.colorbar(image, ax=ax, shrink=.8, orientation ="horizontal")
    else:
        plt.title("Power Spectrum for Pixel %d" % args.pixel)

    #save with filename generated by args
    if args.save:
        name = "plaw-%s-s%d-b%d-l%d-u%d" % (args.file[:-4], args.smoothing, 
               args.bandsize, args.lower, args.upper) 

        if args.moll:
            name += "-moll"

        #save in specified location
        plt.savefig(args.save + name)
        print "Saved as %s.npy in %s" % (name, args.save)

    #or show
    else:
        plt.show()
