#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
import time

desc = "Takes a .npy 2d histogram output by extractor and calculates" \
     +  "the best fit power law for each pixel, plots them"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str,
                    help = ".npy file to be read")
parser.add_argument("-s","--smoothing", type = float, default = 0., 
                    help = "smoothing radius in degrees, default to none")
parser.add_argument("-n", "--normed", default = False, action = "store_true",
                    help = "Whether to renormalize the data after smoothing")
parser.add_argument("-m", "--merge", type = int, default = 0, 
                    help = "number of NChannel bins to merge together")
parser.add_argument("-o", "--outfile", type = str, default = "",
                    help = "name of saved file, default to original name with"+
                    "parameters appended")
args = parser.parse_args()

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
    print "Averaging..."
    zenavg = np.zeros((len(zbands), BINS))
    pavgs = np.zeros(len(zbands))

    #add the pixels in each bin together for sum, calculate best fit for each
    for i, zenbin in enumerate(zbands):
        zcut = (thetas == zenbin)
        zenavg[i] += np.sum(matrix[zcut], axis=0)/float(np.sum(zcut))

    #if true bins, reassign bin width and bins
    if BINS == 1000:
        binwidth = np.zeros(matrix.shape[1]) + 1 
        ybins = np.linspace(0,1000,1000)

    #if true bins but merged, still reassign
    elif args.merge and BINS == 1000/args.merge:
        binwidth = np.zeros(matrix.shape[1]) + 1 
        ybins = np.linspace(0,1000,1000/args.merge)

    print "Saving..."
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = args.file[:-4]

    if args.smoothing:
        outfile += "-s%d" % args.smoothing
    if args.merge:
        outfile += "-m%d" % args.merge
    if args.normed:
        outfile += "-n"
    np.save(outfile, matrix)

    print "Time taken: %d seconds" % (time.time() - start)
