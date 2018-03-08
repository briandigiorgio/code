#!/usr/bin/env python

import numpy as np
import argparse
import glob

desc = "takes a file or directory of .npy files with backwards theta values" +\
       " and corrects them"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str, default = "",
                    help = "file to reverse if only single file")
parser.add_argument("-d", "--directory", type = str, default = "",
                    help = "directory to look for files")
parser.add_argument("-o", "--outfile", type = str, default = "",
                    help = "outfile to write reversed file to, default to" +
                    "writing over old file")
args = parser.parse_args()

#take a matrix and flip it around on the x axis
def reverse(matrix):
    newmatrix = np.zeros(matrix.shape)
    for i in range(matrix.shape[0]):
        newmatrix[i] = matrix[-(i+1)]
    return newmatrix

if __name__ == "__main__":

    #if only for single file, load it and reverse it
    if args.file:
        newmatrix = reverse(np.load(args.file))

        #set appropriate outfile and save it
        if args.outfile:
            outfile = args.outfile
        else:
            outfile = args.file

        np.save(outfile, newmatrix)

    #if whole directory, make list of files in directory
    elif args.directory:
        files = glob.glob(args.directory)
        files.sort()

        #for each file, load it and reverse it
        counter = 0
        for file in files:
            newmatrix = reverse(np.load(file))

            #give it the appropriate name and save it
            if args.outfile:
                outfile = args.outfile + str(counter)
            else:
                outfile = file
            np.save(outfile, newmatrix)
            counter += 1

    else:
        print "Please specify either a file or a directory to reverse"
