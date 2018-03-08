#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import argparse

desc = "Takes a .npy 2d histogram output by extractor.py and plots it"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str, help = ".npy file to be read")
parser.add_argument("-c", "--compression", type = int, default = 512, help = "how much the pixels should be compressed in the graph, should be a power of 2 (default 512)")
parser.add_argument("-a", "--alternate", type = bool, default = False, help = "a parameter than can fix the program if the 'needs more than one value to unpack' error comes up")
args = parser.parse_args()

if __name__ == "__main__":
    #read the matrix in from the .npy file, define parameters
    #has two versions because of some error I can't figure out
    if args.alternate:
        matrix = np.load(args.file)[0]
    else:
        matrix = np.load(args.file)

    npix = matrix.shape[0]
    c = args.compression
    xrange = npix/c
    
    #define new matrix with rows of matrix collapsed down into bins
    b = np.array([np.sum(matrix[c*i:c*(i+1)], axis=0) for i in range(xrange)])

    #plot new matrix 
    fig, ax = plt.subplots()
    cax = ax.matshow(np.log10(b.T), origin = "lower")

    plt.xlabel("Pixels (Bin size: ~%d)" % (49152/xrange))
    plt.ylabel("NChannel (exponential bins)")
    ax.set_title("NChannel by Pixels with Log Counts")

    cbar = fig.colorbar(cax)
    ax.xaxis.tick_bottom()
    plt.show()

