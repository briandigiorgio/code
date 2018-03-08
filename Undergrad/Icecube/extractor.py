#!/usr/bin/env python

import ROOT, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import healpy as hp
import time
from tempfile import TemporaryFile
import argparse
import resource

def using(point = ""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return ("MB of memory used: " + str((usage[2]*resource.getpagesize()/1000000.)))

desc = "Imports root files from a directory and extracts all of the events with their absolute coordinates and nchannel, histograms them and saves the resulting matrix. Also plots it"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-d", "--directory",  type = str, 
       help = 'specify directory to look for root files')
parser.add_argument("-f", "--filenum", type = int, default = 0, 
       help = "index of the file that you want to look at in the directory")
parser.add_argument("-o","--outfile", type = str, default = "outfile",
       help = "file that the matrix is saved to")
parser.add_argument("-y", "--yaxis", type=str, choices=["exp","rawexp","true"],
       help = "scale of y axis (exp, rawexp, or true)")
parser.add_argument("-s", "--start", type = int, default = 0,
       help = "event number to start at, default 0")
parser.add_argument("-e", "--end", type = int, default = 0,
       help = "event number to end, default: end file")
args = parser.parse_args()

if __name__ == "__main__":
    #imports all files from the specified directory with correct format
    files = glob.glob('%s/ic*.root' % args.directory)
    files.sort()

    #defines parameters
    nevents = args.start
    NSIDE = 64
    NPIX = hp.nside2npix(NSIDE)

    #define parameters, cutoff ignores pixels that will be masked
    cutoff = 34688
    npix = [] 
    nchan = []
    nx = NPIX - cutoff
    ny = 100
    xbins = np.linspace(0, nx, nx+1)
    print "###################################################################"
    print
    print "Extracting ROOT file number %d in %s..." \
          % (args.filenum, args.directory)

    ###Generating y bins in different ways depending on user input
    # modified exponential growth, bin sizes grow in discrete steps (default)
    if args.yaxis == "exp":
        print "Generating with modified exponential y axis..."
        yrange = np.arange(ny-1)
        ybins = [6.,7.]
        offset = 0.
        step = 1.

        #generate  50 bins ranging from 0 to 1000 that grow exponentially
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

    #unaltered exponential axis, more true than modified one but less useful
    if args.yaxis == "rawexp":
        print "Generating with raw exponential y axis..." 
        yrange = np.arange(51)
        ybins = [6.,7.]
        offset = 0.
        step = 1.

        for i in yrange:
            newbin = np.trunc(np.exp(np.log(1000) * i/50)) + offset
            ybins.append(newbin)
            diff = ybins[-1] - ybins[-2]

            while diff < 1:
                ybins[-1] += 1
                offset += 1
                diff += 1

    #true y axis, equal size bins up to 1000
    if args.yaxis == "true":
        print "Generating with true y axis..."
        ybins = np.linspace(0,1000,1001)

    #make the last bin an overflow bin
    ybins[-1] = 10000
    print 

    #create a matrix to save that doesn't include masked pixels
    savematrix = np.zeros((nx, len(ybins)-1))

    #opens specified file, gets nchannel, zenith, and az angle for each event
    #converts to healpy pixel, adds all to arrays
    start = time.time()
    f = ROOT.TFile(files[args.filenum])
    t = f.CutDST
    nentries = t.GetEntriesFast()

    #set end of the range if specified
    if not args.end:
        nrange = nentries
    else:
        nrange = args.end

    #helpful print statements
    print "There are %d events in" % nentries, files[args.filenum]
    if args.start or args.end:
        print "But you are running from %d to %d" % (args.start, args.end)
    print 
    print "Looping..."

    #use appropriate while loop to loop over desired amount of events
    while nevents < nrange:

        #extract the data from the entry
        t.GetEntry(nevents)
        theta = (90- t.DecDeg) * np.pi/180

        #if the event not in the masked pixel, cache its data in a list
        if theta > 2.:
            phi = t.RADeg * np.pi/180
            npix += [hp.ang2pix(NSIDE, theta, phi)-cutoff]
            nchan += [t.NChannels]

        nevents += 1

        #clear cache by writing the stored values to a histogram
        if nevents % 100000 == 0:
            savematrix += np.histogram2d(npix,nchan,bins=[xbins,ybins])[0]

            if nevents % 1000000 == 0:
                print "Read %d million events, " % (nevents/1000000), using()

            #clear cache
            npix = []
            nchan = []


    end = time.time()
    print "Time for %d events: %f seconds" % (nevents, end-start)
    print

    #save a 2D histogram to the file specified by user
    print "Saving..."
    np.save(args.outfile, savematrix)
    print using()

    print "Done."
    print 
    print "###################################################################"
