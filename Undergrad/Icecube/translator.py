#!/usr/bin/env python

import numpy as np
import argparse
import glob
import healpy as hp

desc = "take a output histogram from extractor that has been truncated and" \
       "and expand it back out to 49152"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str, default = "",
                    help = "file to translate")
parser.add_argument("-d", "--directory", type = str, default = "",
                    help = "whole directory to translate, overriden by file")
parser.add_argument("-o", "--outfile", type = str, default = "",
                    help = "outfile for translated matrix, default to writing" 
                    + "over old file")
args = parser.parse_args()

#take a matrix that has been truncated and restore it to 49152 width
def translate(matrix):

    #define new matrix with correct dimensions, find difference in width
    width = matrix.shape[0]
    newmatrix = np.zeros((49152, matrix.shape[1])) + hp.UNSEEN
    shift = 49152 - width
    
    #populate the matrix in the correct cells
    for i in range(width):
        newmatrix[i + shift] = matrix[i]
    
    return newmatrix

if __name__ == "__main__":

    #if only one file, load it and translate it
    if args.file:
        newmatrix = translate(np.load(args.file))

        #save the file with correct outfile name
        if args.outfile:
            outfile = args.outfile
        else:
            outfile = args.file

        np.save(outfile, newmatrix)

    #if whole directory, load a list of all of the files
    elif args.directory:
        files = glob.glob(args.directory)
        files.sort()

        #translate each file and save it with correct name
        for file in files:
            newmatrix = translate(np.load(file))

            if args.outfile:
                outfile = args.outfile
            else:
                outfile = file 

            np.save(outfile, newmatrix)
