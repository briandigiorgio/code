#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
import time
import matplotlib.pyplot as plt
from scipy import optimize as sp
from scipy import stats

desc = "Takes a .npy 2d histogram output by extractor and calculates" \
     + "the the two factor Kolmogorov-Smirnov statistic for each pixel and" \
     + "the average of its declination band, plots results"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str,
                    help = ".npy file to be read")
parser.add_argument("-s","--smoothing", type = float, default = 0., 
                    help = "smoothing radius in degrees, default to none")
parser.add_argument("-n", "--normed", default = False, action = "store_true",
                    help = "Whether to renormalize the data after smoothing")
parser.add_argument("-b", "--bandsize", type = int, default = 0, 
                    help = "the size of the declination bands in pixels, "
                    + "odd numbers make more sense")
parser.add_argument("-m", "--merge", type = int, default = 0, 
                    help = "number of NChannel bins to merge together")
parser.add_argument("--min", type = float, help = "minimum plot value")
parser.add_argument("--max", type = float, help = "maximum plot value")
parser.add_argument("--moll", default = False, action = "store_true", 
                    help = "plot in mollview instead of orthview")
parser.add_argument("-k", "--ks", default = False, action = "store_true",
                    help = "show a plot of ks values instead of p values")
parser.add_argument("-l", "--lower", type = int, default = 0, 
                    help = "lower cut of bins to find fitline of")
parser.add_argument("-u", "--upper", type = int, default = -1,
                    help = "upper cut of bins to find fitline of")
parser.add_argument("--save", type = str, default = "", 
                    help = "saves plot to filepath given instead of showing it")
parser.add_argument("--mask", type = int, default = 0,
                    help = "number of degrees to mask around the pole")
args = parser.parse_args()

def kstest(d1,d2):
    d1sum = np.sum(d1)
    d2sum = np.sum(d2)
    length = len(d1)
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
   
    #merge NChannel bins together according to how big the user wants the bins
    if args.merge:
        print "Merging bins..."
        mmatrix = np.zeros((NPIX, BINS/args.merge))

        #divide the index of each bin by merging factor and put into new bin
        for pixel in range(NPIX):
            for bin in range(BINS):
                mindex = bin/args.merge
                if mindex >= BINS/args.merge:
                    mindex -= 1

                mmatrix[pixel][mindex] += matrix[pixel][bin]

        #reassign variables
        matrix = mmatrix
        BINS = BINS/args.merge

    #smooths the data over some angular radius that user specifies
    #does nothing if no radius is given
    if args.smoothing:
        print "Smoothing..."
        smatrix = np.zeros(matrix.shape)

        #convert to radians
        smoothrad = args.smoothing * np.pi/180

        #find the neighbors of the pixel, add values together to make new value
        for pixel in range(NPIX):
            vec = hp.pix2vec(NSIDE, pixel)
            neighbors = hp.query_disc(NSIDE, vec, smoothrad)
            smatrix[pixel] += matrix[neighbors].sum(axis = 0)

        #renormalize the values by dividing by number of pixels combined
            if args.normed:
                smatrix[pixel] = smatrix[pixel]/float(len(neighbors))

        #reassign variables
        matrix = smatrix

    #make arrays of sums for each bin and the power law exponent in each band
    print "Extracting..."
    zenavg = np.zeros((len(zbands), BINS))
    pavgs = np.zeros(len(zbands))

    #add the pixels in each bin together for sum
    for i, zenbin in enumerate(zbands):
        zcut = (thetas == zenbin)
        zenavg[i] += np.sum(matrix[zcut], axis=0)/float(np.sum(zcut))

    #make arrays for storing results
    ksvals = np.zeros(NPIX)
    pvals = np.zeros(NPIX)
    maxbin = np.zeros(NPIX)

    #if true bins, reassign bin width and bins
    if BINS == 1000:
        binwidth = np.zeros(matrix.shape[1]) + 1 
        ybins = np.linspace(0,1000,1000)

    #if true bins but merged, still reassign
    elif args.merge and BINS == 1000/args.merge:
        binwidth = np.zeros(matrix.shape[1]) + 1 
        ybins = np.linspace(0,1000,1000/args.merge)

    lower,upper=np.digitize([args.lower,args.upper],ybins)

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

                average = average / args.bandsize

            #if value hasn't been calculated already, run through same process
            elif not pavgs[zenbin]:
                average = np.zeros(len(zenavg[zenbin]))
                average += zenavg[zenbin]

            #finds difference in power law exponent for each pixel
            if args.lower or args.upper != -1:
                pixcut = matrix[pixel][lower:upper]
                avgcut = average[lower:upper]
                #ksvals[pixel], pvals[pixel] = stats.ks_2samp(pixcut,avgcut)
                ksvals[pixel],pvals[pixel],maxbin[pixel]=kstest(pixcut,avgcut)
            else:
                #ksvals[pixel], pvals[pixel]=stats.ks_2samp(np.log10(
#matrix[pixel]),np.log10(average))
                ksvals[pixel],pvals[pixel],maxbin[pixel]=kstest(matrix[pixel],
                    average)

            if not pixel % 1000 and pixel:
                print "Pixels processed: %d,000" % (pixel/1000)
                print ksvals[pixel], pvals[pixel], maxbin[pixel]

    print "Time taken: %f seconds" % (time.time() - start)

    #hide pixels in the northern hemisphere and plot 
    print "Plotting..."

    #hide pixels that are more than 65 deg from pole or are too close to pole
    for i in range(NPIX):
        if thetas[i] < 2:
            ksvals[i] = hp.UNSEEN
            pvals[i] = hp.UNSEEN

        #mask weirdness at pole
        if args.mask and (np.pi - thetas[i]) * 180/np.pi < args.mask:
            ksvals[i] = hp.UNSEEN
            pvals[i] = hp.UNSEEN

    #show in IC standard rotated half sky mollview
    if args.moll and args.ks:
        hp.mollview(maxbin,min=args.min,max=args.max,cbar=False,rot=180)
        plt.ylim(-1, .005)

    elif args.moll:
        hp.mollview(pvals,min=args.min,max=args.max,cbar=False,rot=180)
        plt.ylim(-1, .005)

    #show in normal IceCube orthview of southern hemisphere if not specified
    elif args.ks:
        hp.orthview(maxbin, rot = [0,-90,180], half_sky = True, min = args.min,
                    max = args.max, cbar = False)
    else:
        hp.orthview(pvals, rot = [0,-90,180], half_sky = True, min = args.min,
                    max = args.max, cbar = False)

    #gridlines
    hp.graticule()

    #a whole long thing to generate a caption for the plot based off args
    print
    title = "Two Factor Kolmogorov-Smirnov plot for %s showing" % args.file 

    if args.ks:
        title += " KS values"
    else:
        title += " p values"
        
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
    if args.normed:
        title += ", normalized"

    #make each line of the caption the correct length
    for i in range(len(title)/80):
        title = title[:79 * (i+1)] + "- \n" + title[79 * (i+1):]

    #add title and caption
    print title
    #plt.figtext(.1, .02, title)
    plt.title("Kolmogorov-Smirnov Test")
     
    #make color bar correct size
    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]
    cbar = fig.colorbar(image, ax=ax, shrink=.8, orientation ="horizontal")
 
    #save with filename generated by args
    if args.save:
        if args.ks:
            type = "ks"
        else:
            type = "pval"

        name = "ks-%s-%s-s%d-b%d-l%d-u%d" % (args.file[:-4], type, 
               args.smoothing, args.bandsize, args.lower, args.upper) 

        if args.moll:
            name += "-moll"

        #save in specified location
        plt.savefig(args.save + name)
        print "Saved as %s.npy in %s" % (name, args.save)

    #or show
    else:
        plt.show()
