#!/usr/bin/env python

import os
import time

'''
This is kind of an ad hoc program that I wrote to generate a whole bunch of
of plots. It's not really supposed to be usable with arguments and stuff
because typing out that so many arguments that are so long would get super
annoying very quickly. But anyways, if you just put whatever you want in
thefiles, smoothing levels, bandsizes, and limits in the things below, it'll
run bayes.py and powerlaw.py on all of them. It takes a long time.
'''

if __name__ == "__main__":
    start = time.time()

    #parameters to run over, it will run all possible combinations
    #choose carefully, it snowballs very quickly
    files = ["revall.npy","revIC86.npy","truetIC86.npy","truetall.npy"]
    smoothing = [0,5,20]
    bandsize = [0,5,11]
    limits = [(33,100)]

    #where to save the plots
    save = "~/plots/"

    command = ""
    for f in files:
        for s in smoothing:
            for b in bandsize:
                for l in limits:

                    #generate the actual terminal commands
                    pcommand = "./powerlaw.py -f %s -s %d -b %d -l %d -u %d " \
                             % (f, s, b, l[0], l[1]) \
                             + "--min -.01 --max .01 --mask 12 --save " + save 

                    bcommand = "./bayes.py -f %s -s %d -b %d -l %d -u %d " \
                             % (f, s, b, l[0], l[1]) \
                             + "--mask 12 --save " + save

                    #run each command in orthview and mollview
                    #run each command with orthview and mollview
                    print pcommand
                    os.system(pcommand)
                    print
                    print pcommand + " --moll"
                    os.system(pcommand + " --moll")
                    print

                    print bcommand
                    os.system(bcommand)
                    print
                    print bcommand + " --moll"
                    os.system(bcommand + " --moll")
                    print

    print "Total time taken: ", time.time() - start
