#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
from scipy import special as sp
import time
import matplotlib.pyplot as plt

desc = "Takes a .npy 2d histogram output by extractory and calculates" \
     +  "the Bayes factor for each pixel compared to its dec band, plots them"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str,
                help = ".npy file to be read")
parser.add_argument("-a", "--alternate", default = False, action = "store_true",
               help = "a parameter than can fix the program if the " +
               "'needs more than one value to unpack' error comes up")
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
parser.add_argument("-l", "--lower", type = int, default = 0, 
                    help = "lower cut of bins to find fitline of")
parser.add_argument("-u", "--upper", type = int, default = -1,
                    help = "upper cut of bins to find fitline of")
parser.add_argument("--mask", type = int, default = 0,
                    help = "number of degrees to mask around the pole")
parser.add_argument("--save", type = str, default = "", 
                    help = "saves plot to filepath given instead of showing it")
args = parser.parse_args()

if __name__ == "__main__":
    start = time.time()

    print "Importing..."
    #dumb if statement to take care of errors
    if args.alternate:
        matrix = np.load(args.file)[0]
    else:
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
        if theta != zbands[band]:
            band += 1
        binnum[pixel] = band

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

    #make arrays of sums for each bin and the total sum in each band
    print "Extracting..."
    if args.upper == -1:
        args.upper = len(ybins) - 1
    lower,upper=np.digitize([args.lower,args.upper],ybins)

    #add the pixels in each bin together for sum, add across all bins for total
    zentots = np.zeros(len(zbands))
    zensum = np.zeros((len(zbands), upper-lower))
    for i, zenbin in enumerate(zbands):
        zcut = (thetas == zenbin)
        zensum[i] += np.sum(matrix[zcut], axis=0)[lower:upper]/np.sum(zcut)
    zentots = np.sum(zensum, axis=1)

    #apply formula from paper to find lnB for each pixel
    print "Calculating..."
    lnB = np.zeros(NPIX)
    for pixel,zenbin in enumerate(binnum):

        #only look at pixels that are less than 65 deg from pole
        if thetas[pixel] > 2:
            gammasum = np.zeros(NPIX)
            average = np.zeros(len(zensum[zenbin]))
            average += zensum[zenbin]
            total = 0
            total += zentots[zenbin]

            #widens the declination bands used in averages in the formula
            if args.bandsize:
                size = 1
                delta = 1

                #while the bands are not as big as you want them, widen them
                #adds a row of pixels above then below then above then below...
                while size < args.bandsize:

                    #reverses band expansion at pole so it doesn't go over
                    if zenbin + delta >= 255:
                        delta = -delta - 1

                    #actually adds the rows to the necessary values
                    average += zensum[zenbin + delta]
                    total += zentots[zenbin + delta]

                    #adjusts offset to look in the correct place for next band
                    size += 1
                    delta = -delta
                    if delta > 0:
                        delta += 1 

                #renormalize counts
                average = average/args.bandsize
                total = total/args.bandsize

            #formula from paper applied to each pixel (arxiv.org/abs/1106.1392)
            f1 = matrix[pixel][lower:upper]
            gammasum[pixel] = np.sum(sp.gammaln(f1 + 1) \
                            + sp.gammaln(average + 1) \
                            - sp.gammaln(f1 + average + 2))

            lnB[pixel] = gammasum[pixel] \
                       +  sp.gammaln(np.sum(f1) + total + 2) \
                       - sp.gammaln(np.sum(f1) + 1) \
                       - sp.gammaln(total + 1)

            if pixel % 1000 == 0:
                print "Pixels processed: %s,000" % str(pixel/1000)

        #change into log10
        logB = lnB/np.log(10)

    print "Time taken: %f seconds" % (time.time() - start)

    #hide pixels more than 65 deg from pole and also right next to pole
    print "Plotting..."
    for i in range(NPIX):
        if thetas[i] < 2:
            logB[i] = hp.UNSEEN

        if args.mask and (np.pi - thetas[i]) * 180/np.pi < args.mask:
            logB[i] = hp.UNSEEN

    #show plot in standard IC rotated mollview half plot
    if args.moll:
        hp.mollview(logB, min = args.min, max = args.max, rot = 180,
                    cbar = False)
        plt.ylim(-1, .005)

    #show plot in standard IC rotated orthview with only southern sky
    else:
        hp.orthview(logB, rot = [0,-90,180], half_sky = True, min = args.min,
                    max = args.max, cbar = False)
    #standard grid lines
    hp.graticule()

    #whole complicated thing to generate a caption for the plot based on args
    title = "Bayes Factor plot for " + args.file
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

    #make the caption have 80 characters per line
    for i in range(len(title)/80):
        title = title[:80 * (i+1) - 1] + "- \n" + title[80 * (i+1) - 1:]

    print title

    #put title and caption on plot
    plt.title("Bayes Factor")
    plt.figtext(.1, .02, title)

    #get color bar correct size, place, and ticks
    fig = plt.gcf()
    ax = plt.gca()
    image = ax.get_images()[0]
    cbar = fig.colorbar(image, ax = ax, shrink = .8, orientation = "horizontal")

    #save plot if requested with title generated from args
    if args.save:
        name = "bayes-%s-s%d-b%d-l%d-u%d" % (args.file[:-4], args.smoothing, 
               args.bandsize, args.lower, args.upper) 

        if args.moll:
            name += "-moll"

        #save in directory specified
        plt.savefig(args.save + name)
        print "Saved as %s.npy in %s" % (name, args.save)

    #or just show it
    else:
        plt.show()
