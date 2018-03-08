#!/usr/bin/env python

import numpy as np
import argparse
import glob
import time 
import resource

def using(point = ""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return ("MB of memory used: " + str((usage[2]*resource.getpagesize()/1000000.)))

desc = "Combines the .npy histograms output by extractor in a directory into" \
     + "one grand histogram"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-d", "--directory", type = str,
                    help = "directory of files to combine")
parser.add_argument("-n", "--number", type = int, default = 0,
                    help = "number of files to combine. Defaults to all")
parser.add_argument("-o", "--outfile", type = str, default = "combined",
                    help = "name of the combined .npy file")
parser.add_argument("-a", "--alternate", default = False, action ="store_true",
                    help = "set to True if the 'needs more than one value to " \
                    + "unpack' error comes up")
parser.add_argument("-s", "--string", type = str, default = "extracted",
                    help = "beginning of filename to draw all files from")
args = parser.parse_args()

if __name__ == "__main__":
    #load all files in directory with specific filenames
    start = time.time()
    print "Importing..."
    files = glob.glob(args.directory + args.string + "*.npy")
    files.sort()

    #load specified number of files    
    if not args.number:
        numfiles = len(files)
    else:
        numfiles = args.number

    #load the first file
    if args.alternate:
        matrix = np.load(files[0])[0]
    else:
        matrix = np.load(files[0])
    
    print "Reading..."
    print "1 file read...", using()
    n = 1

    #load the rest of the files and add them to the first one
    while n < numfiles:
        file = np.load(files[n])
        matrix += file
        n += 1
        print "%d files read..." % n, using()

    #save the matrix
    np.save(args.directory + args.outfile, matrix)
    print "There are now %d events in the file" % matrix.sum()
    print "Done. Time taken: %f," % (time.time() - start), using()
