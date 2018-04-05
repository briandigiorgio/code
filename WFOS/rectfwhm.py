#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from optexttools import *
import time

if __name__ == "__main__":
    start = time.time()
    npix = 180
    size = 9/npix
    pixels = 1/size
    ff = 10000
    fwhms = np.linspace(.3,.9,10)
    noises = [1000]
    shape = (len(noises), len(fwhms))
    cmap = make_cmap(len(noises), 'inferno')
    nruns = 50

    flux = np.zeros(shape)
    yerr = np.zeros(shape)
    sn = np.zeros(shape)
    aflux = np.zeros(shape)
    ayerr = np.zeros(shape)
    asn = np.zeros(shape)

    sns = np.zeros(nruns)
    asns = np.zeros(nruns)
    errors = np.zeros(nruns)
    aerrors = np.zeros(nruns)

    counter = 0
    total = len(noises)*len(fwhms)
    progtime = time.time()
    progressbar(counter, total)

    for i in range(len(noises)):
        for j in range(len(fwhms)):
            noise = noises[i]/np.sqrt(pixels)
            for k in range(nruns):
                flux[i][j], errors[k], sns[k] = optrext(size, fwhm = fwhms[j],
                                        npix = npix, ff = ff, noise = noise)
                aflux[i][j], aerrors[k], asns[k] = optrext(size, fwhm=fwhms[j],
                            npix = npix, ff = ff, noise = noise, a = True)

            sn[i][j] = np.average(sns)
            asn[i][j] = np.average(asns)
            #yerr[i][j] = np.sqrt(np.sum(np.power(errors/nruns, 2)))
            #ayerr[i][j] = np.sqrt(np.sum(np.power(aerrors/nruns, 2)))

            counter += 1
            progtime = progressbar(counter, total, progtime)

        plt.plot(fwhms, sn[i], c = cmap[i], label = 'Optimal Extraction')
        plt.plot(fwhms, asn[i], c = cmap[i], ls = '--',
            label = 'Aperture Sum')
            #plt.fill_between(fibers, sn[i]-yerr[i], sn[i]+yerr[i], facecolor =
            #        cmap[i], alpha = .2)

    print(asn/sn)
    plt.title("S/N Ratios for Varying FWHM with Source/Background = 100")
    plt.xlabel("PSF FWHM (arcsec)")
    plt.ylabel("S/N")
    plt.grid(True)
    plt.legend(title = "Background/Source")
    print("\nTime taken: ", time.time()-start)
    plt.show()
