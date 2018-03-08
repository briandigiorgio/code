#!/usr/bin/env python

import numpy as np
import healpy as hp
import argparse
import matplotlib.pyplot as plt

desc = "plot the counts per pixel of an extractor output .npy file"
parser = argparse.ArgumentParser(description = desc)
parser.add_argument("-f", "--file", type = str, help = "file to read")
parser.add_argument("--moll", default = False, action = "store_true")
parser.add_argument("--orth", default = False, action = "store_true")
parser.add_argument("-l", "--log", default = False, action = "store_true")
args = parser.parse_args()

if __name__ == "__main__":
    #load file, sum all counts for each pixel
    matrix = np.load(args.file)
    sums = np.sum(matrix, axis = 1)

    if args.log:
        sums = np.log10(sums)

    #reset unseen values
    for i in range(len(sums)):
        if sums[i] < 0:
            sums[i] = hp.UNSEEN

    #plot it how the user wanted it to
    if args.moll:
        hp.mollview(sums)

    elif args.orth:
        hp.orthview(sums, rot  = [0,-90,0])

    else:
        for i in range(34688):
            sums[i] = hp.UNSEEN

        hp.orthview(sums, rot = [0,-90,180], half_sky = True)

    plt.show()
