#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from optexttools import *
import time

if __name__ == "__main__":
    start = time.time()
    fibers = [1,2,4,8,16,32]
    ff = 10000
    #ff = [100,1000,10000,100000]
    noises = [1,100,1000,10000,100000]
    #noises = 1000
    shape = (len(noises), len(fibers))
    #shape = (len(ff), len(fibers))
    cmap = make_cmap(len(noises), 'inferno')
    #cmap = make_cmap(len(ff), 'inferno')
    nruns = 100
    normed = False
    fwhm=.8

    flux = np.zeros(shape)
    yerr = np.zeros(shape)
    sn = np.zeros(shape)
    aflux = np.zeros(shape)
    ayerr = np.zeros(shape)
    asn = np.zeros(shape)
    snnorm = np.zeros(shape)
    asnnorm = np.zeros(shape)
    sns = np.zeros(nruns)
    asns = np.zeros(nruns)
    errors = np.zeros(nruns)
    aerrors = np.zeros(nruns)

    counter = 0
    total = len(noises)*len(fibers)
    progtime = time.time()

    for i in range(len(noises)):
    #for i in range(len(ff)):
        for j in range(len(fibers)):
            noise = noises[i]/(fibers[j])
            #noise = noises/(fibers[j]**2)
            for k in range(nruns):
                flux[i][j], errors[k], sns[k] = optext(fibers[j], fwhm=fwhm,
                                                 ff = ff, noise = noise)
                aflux[i][j], aerrors[k], asns[k] = optext(fibers[j], fwhm=fwhm,
                            ff = ff, noise = noise, a = True)

            sn[i][j] = np.average(sns)
            asn[i][j] = np.average(asns)
            snnorm[i][j] = sn[i][j]*np.sqrt(ff+noises[i]**2)/ff
            asnnorm[i][j] = asn[i][j]*np.sqrt(ff+noises[i]**2)/ff
            yerr[i][j] = np.sqrt(np.sum(np.power(errors/nruns, 2)))
            ayerr[i][j] = np.sqrt(np.sum(np.power(aerrors/nruns, 2)))

            counter += 1
            #print("Done %d/%d" % (counter, total))
            progtime = progressbar(counter, total, progtime)
            
        if normed:
            plt.loglog(fibers, snnorm[i], c = cmap[i], label = noises[i])
            plt.loglog(fibers, asnnorm[i], c = cmap[i], ls = '--')
        else:
            plt.loglog(fibers, sn[i], c = cmap[i], label = (noises[i]**2/ff))
            plt.loglog(fibers, asn[i], c = cmap[i], ls = '--')
            #plt.fill_between(fibers, sn[i]-yerr[i], sn[i]+yerr[i], facecolor =
            #        cmap[i], alpha = .2)

    plt.title("S/N Ratios for Different Background Noise for FWHM = " +
            str(fwhm) + '"')
    plt.xlabel("Number of Fibers on a Side")
    plt.ylabel("S/N")
    plt.grid(True)
    plt.legend(title = "Background/Source")
    #plt.savefig("fluxcomp.png")
    print("\nTime taken: ", time.time()-start)
    #print(yerr)
    plt.show()
