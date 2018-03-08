#!/usr/bin/env python

import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stats
from scipy.signal import convolve2d
from matplotlib import colors, rc, cm
from astropy.modeling import functional_models as dists

def gaussian(x, mean, std):
    return (1/(std*np.sqrt(2*np.pi)))**np.exp(-1*((x-mean)**2)/(2*std**2))

def addnoise(psf, noise):
    return psf + np.random.normal(0, noise, size = psf.shape)

#returns a moffat psf centered at [x,y] = center evaluated at pos (2darray)
#beta = 2.9 is given by Gemini GLAO in 
#http://www.ctio.noao.edu/~atokovin/papers/Gemini-GLAO.pdf
#parameter naming is messed up because astropy uses unconventional names
def moffat(center, x, y, fwhm, beta = 2.9):
    alpha = fwhm/(2*np.sqrt(2**(1/beta)-1))
    tpsf = np.zeros((len(x),len(y)))
    for i, xi in enumerate(x):
        for j, yj in enumerate(y):
            tpsf[i,j] = dists.Moffat2D.evaluate(xi, yj, 1, x_0 = center[0],
                    y_0 = center[1], gamma = alpha, alpha = beta)
    return tpsf

def optimalextract(data, tdata, var, aperture = False):
    #method from Naylor 1998 Eqs. 9-12
    #data = observed counts in fibers, tdata = theoretical psf values
    #returns flux, error on flux measurement, S/N ratio
    sumvar = np.sum(np.square(tdata)/var)
    weights = (tdata/var)/sumvar

    if aperture:
        weights = np.ones((data.shape))

    flux = np.sum(weights * data)
    error = np.sqrt(np.sum(np.square(weights) * var))
    #sn = flux * np.sqrt(sumvar)
    sn = flux / error
    return flux, error, sn 

def apsum(data, var):
    flux = np.sum(data) 
    error = np.sqrt(np.sum(var))
    sn = flux/error
    return flux, error, sn

def optext(nfibers, fwhm, npix = 128, noise = 50, ff = 10000, a = False, m =
        True):
    #define variables
    std = fwhm/2.355 #convert fwhm to stdev
    fibersize = npix//nfibers #size of fiber in pixels
    fiber = np.ones((fibersize,fibersize)) #square mask to represent pixel
    size = 1

    #make grid of positions across pixel space
    x = np.linspace(-.5*size,.5*size,npix)
    y = np.linspace(-.5*size,.5*size,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    #define 2d pdf, evaluate to make a grid of values, normalize
    if m:
        tpsf = moffat((0,0), x, y, fwhm, 2.9)
        tpsf = tpsf/np.sum(tpsf)
    else:
        pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
        tpsf = pdf.pdf(pos)

    tpsf = tpsf/np.sum(tpsf)
    psf = tpsf * ff

    #convolve fibers with psfs to get simulated fiber data
    alldata = convolve2d(psf, fiber, mode = 'valid')
    talldata = convolve2d(tpsf, fiber, mode = 'valid')

    #select only the data that would correspond to unoverlapping fibers
    xpos = list(range(0,alldata.shape[0]-fibersize,fibersize))\
        +[alldata.shape[0]-1]
    ypos = list(range(0,alldata.shape[1]-fibersize,fibersize))\
        +[alldata.shape[1]-1]
    data = np.random.poisson(alldata[xpos,:][:,ypos])
    tdata = talldata[xpos,:][:,ypos]
    data = addnoise(data, noise)
    var = np.square(noise) + data #variance of noise
    '''
    plt.imshow(data, cmap = "inferno")
    plt.title("Noise = %d" % noise)
    plt.show()
    '''
    #if a:
    #    return apsum(data, var)
    return optimalextract(data, tdata, var, aperture = a) 
    #return var


def difftest(nfibers, fwhm, npix = 128, noise = 50, ff = 10000, a = False):
    #define variables
    std = fwhm/2.355 #convert fwhm to stdev
    fibersize = npix//nfibers #size of fiber in pixels
    fiber = np.ones((fibersize,fibersize)) #square mask to represent pixel

    #make grid of positions across pixel space
    x = np.linspace(-1,1,npix)
    y = np.linspace(-1,1,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    #define 2d gaussian pdf, evaluate to make a grid of values, normalize
    pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
    tpsf = pdf.pdf(pos)
    tpsf = tpsf/np.sum(tpsf)

    #add noise so it's not perfect
    #psf = ff * addnoise(tpsf, noise)
    psf = tpsf * ff

    #convolve fibers with psfs to get simulated fiber data
    alldata = convolve2d(psf, fiber, mode = 'valid')
    talldata = convolve2d(tpsf, fiber, mode = 'valid')

    #select only the data that would correspond to unoverlapping fibers
    xpos = list(range(0,alldata.shape[0]-fibersize,fibersize))\
        +[alldata.shape[0]-1]
    ypos = list(range(0,alldata.shape[1]-fibersize,fibersize))\
        +[alldata.shape[1]-1]
    data = np.random.poisson(alldata[xpos,:][:,ypos])
    tdata = talldata[xpos,:][:,ypos]
    data = addnoise(data, noise)
    var = np.square(noise) + data #variance of noise

    #plt.imshow(data, cmap = "inferno")
    #plt.colorbar()
    #plt.show()

    flux, error, sn = optimalextract(data, tdata, var, aperture = False)
    aflux, aerror, asn = optimalextract(data, tdata, var, aperture = True)

    return diff(flux, aflux), diff(error, aerror), diff(sn, asn)

def diff(orig, new):
    return ((new - orig)/orig) * 100

#makes an array of length numlines going through the specified color map
#plt.plot(..., c = colors[i]) when plotting multiple lines
def make_cmap(numlines, cmap):
    cnorm = colors.Normalize(vmin = 0, vmax = numlines)
    scalarmap = cm.ScalarMappable(norm = cnorm, cmap = cmap)
    return scalarmap.to_rgba(range(numlines))

def progressbar(progress, total, progtime = 0, tlast = 0):
    end = ''
    print(' '*73, end = '\r')
    if not progtime:
        end = '\r'
    print('Progress: [%s%s] %d%%' % ('▒' * progress, '_' * (total-progress), 
         (progress * 100)//total), end = end)
    if progtime:
        timediff = float(time.time()-progtime)

        if tlast:
            new = np.array(len(tlast))
            for i in range(len(tlast)):
                new[i] = tlast[i]/(total-progress + i)
            new = new.append(timediff)
            print("  Time Remaining: ~%d seconds" %
                   ((total-progress)*np.average(new)), end = '\r')
            return time.time(), new

        print("  Time Remaining: ~%d seconds" %
                ((total-progress)*timediff), end = '\r')
            
        return time.time()

#returns a hexagonal mask array with desired radius
def makehex(rad):
    r3 = np.sqrt(3)
    shape = (2*rad,int(r3*rad))
    hexagon = np.ones(shape)
    for x in range(shape[0]):
        for y in range(shape[1]):
            if y > r3/2*rad + r3*x or y < r3/2*rad - r3*x \
            or y > -r3*x + 5*r3/2*rad or y < r3*x - 3*r3/2*rad:
                hexagon[x,y] = 0
    #print(hexagon)
    #plt.imshow(hexagon)
    #plt.show()
    return hexagon

#convolve a fiber onto a psf at only one specified location
#input center of the fiber as (y,x) coordinate
def convolve(center, fiber, psf):
    odd = (fiber.shape[0]%2, fiber.shape[1]%2) #accounts for odd dimensions
    xoffset = fiber.shape[1]//2
    yoffset = fiber.shape[0]//2
    #multiplies arrays together
    product = fiber * psf[center[0]-yoffset:center[0]+yoffset+odd[0],
            center[1]-xoffset:center[1]+xoffset+odd[1]]
    return np.sum(product)

#performs optimal extraction simulation using a first order hexagonal IFU
#first order meaning 1 outer ring of fibers (7 fiber)
def opthext1(fwhm, npix = 128, size = 1, noise = 1000, ff = 10000, a = False, m
        = True):
    #define variables
    r3 = np.sqrt(3)
    std = fwhm/2.355 #convert fwhm to stdev
    r = int(npix/(3*r3))
    fiber = makehex(r) #square mask to represent pixel
    center = npix//2

    #make grid of positions across pixel space
    x = np.linspace(-.5*size,.5*size,npix)
    y = np.linspace(-.5*size,.5*size,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    if m:
        tpsf = moffat((0,0), x, y, fwhm, 2.9)
        tpsf = tpsf/np.sum(tpsf)
    else:
        pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
        tpsf = pdf.pdf(pos)

    tpsf = tpsf/np.sum(tpsf)
    psf = tpsf * ff
    

    #manually goes through and finds values for all of the fiber locations
    #there's probably a more elegant way to do this
    #letters start at theta=0 and progress counterclockwise
    f0 = convolve((center,center), fiber, psf)
    f1a = convolve((center, int(center+r*r3)), fiber, psf)
    f1b = convolve((int(center+1.5*r), int(center+(r*r3)/2)), fiber, psf)
    f1c = convolve((int(center+1.5*r), int(center-(r*r3)/2)), fiber, psf)
    f1d = convolve((center, int(center-r*r3)), fiber, psf)
    f1e = convolve((int(center-1.5*r), int(center-(r*r3)/2)), fiber, psf)
    f1f = convolve((int(center-1.5*r), int(center+(r*r3)/2)), fiber, psf)

    #put all data into an array as it would be line by line
    #in case I ever figure out how to plot this
    data =   [f1c,f1b,
            f1d,f0,f1a,
              f1e,f1f]
    tdata = data/np.sum(data)

    #add noise, find variance, do optext
    data = addnoise(np.random.poisson(data), noise)
    var = np.square(noise) + data #variance of noise
    return optimalextract(data, tdata, var, aperture = a) 

#same as opthext1 but for a 2nd order bundle (2 rings, 19 fibers)
def opthext2(fwhm, npix = 128, size = 1, noise = 1000, ff = 10000, a = False, m
        = True):
    #define variables
    r3 = np.sqrt(3)
    std = fwhm/2.355 #convert fwhm to stdev
    r = int(npix/(5*r3))
    fiber = makehex(r) #square mask to represent pixel
    center = npix//2

    #make grid of positions across pixel space
    x = np.linspace(-.5*size,.5*size,npix)
    y = np.linspace(-.5*size,.5*size,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    #define 2d pdf, evaluate to make a grid of values, normalize
    if m:
        tpsf = moffat((0,0), x, y, fwhm, 2.9)
        tpsf = tpsf/np.sum(tpsf)
    else:
        pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
        tpsf = pdf.pdf(pos)

    tpsf = tpsf/np.sum(tpsf)
    psf = tpsf * ff
    
    #make all of the different fibers
    #there must be a better way to do this
    f0 = convolve((center,center), fiber, psf)

    f1a = convolve((center, int(center+r*r3)), fiber, psf)
    f1b = convolve((int(center+1.5*r), int(center+(r*r3)/2)), fiber, psf)
    f1c = convolve((int(center+1.5*r), int(center-(r*r3)/2)), fiber, psf)
    f1d = convolve((center, int(center-r*r3)), fiber, psf)
    f1e = convolve((int(center-1.5*r), int(center-(r*r3)/2)), fiber, psf)
    f1f = convolve((int(center-1.5*r), int(center+(r*r3)/2)), fiber, psf)

    f2a = convolve((center, int(center+2*r*r3)), fiber, psf)
    f2b = convolve((int(center+1.5*r), int(center+1.5*r*r3)), fiber, psf)
    f2c = convolve((center+3*r, int(center+r*r3)), fiber, psf)
    f2d = convolve((center+3*r, center), fiber, psf)
    f2e = convolve((center+3*r, int(center-r*r3)), fiber, psf)
    f2f = convolve((int(center+1.5*r), int(center-1.5*r*r3)), fiber, psf)
    f2g = convolve((center, int(center-2*r*r3)), fiber, psf)
    f2h = convolve((int(center-1.5*r), int(center-1.5*r*r3)), fiber, psf)
    f2i = convolve((center-3*r, int(center-r*r3)), fiber, psf)
    f2j = convolve((center-3*r, center), fiber, psf)
    f2k = convolve((center-3*r, int(center+r*r3)), fiber, psf)
    f2l = convolve((int(center-1.5*r), int(center+1.5*r*r3)), fiber, psf)

    #put all of the data together into an array scanning line by line
    #the array is a graphical representation of where all of the fibers are
    data =  [f2e,f2d,f2c,
           f2f,f1c,f1b,f2b,
         f2g,f1d,f0,f1a,f2a,
           f2h,f1e,f1f,f2l,
             f2i,f2j,f2k]
    tdata = data/np.sum(data)

    data = addnoise(np.random.poisson(data), noise)
    var = np.square(noise) + data #variance of noise
    return optimalextract(data, tdata, var, aperture = a)

#one hexagonal fiber, some stuff had to be rewritten so it would work
def opthext0(fwhm, npix = 128, size = 1, noise = 1000, ff = 10000, a = False, m
        = True):
    #define variables
    r3 = np.sqrt(3)
    std = fwhm/2.355 #convert fwhm to stdev
    r = npix//2
    fiber = makehex(r) #square mask to represent pixel
    center = npix//2

    #make grid of positions across pixel space
    x = np.linspace(-.5*size,.5*size,npix)
    y = np.linspace(-.5,.5,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    if m:
        tpsf = moffat((0,0), x, y, fwhm, 2.9)
        tpsf = tpsf/np.sum(tpsf)
    else:
        pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
        tpsf = pdf.pdf(pos)

    tpsf = tpsf/np.sum(tpsf)
    psf = tpsf * ff
    

    #get the actual data
    data = convolve((center,center), fiber, psf)
    data = np.random.poisson(data) + np.random.normal(0, noise)

    #calculate error, s/n
    #doesn't use optext because it isn't necessary and breaks the function
    error = np.square(noise) + data #variance of noise
    sn = data/error
    return data, error, sn

#umbrella function for other opthexts, calls appropriate one
def opthext(fibers,fwhm,npix = 128,size = 1,noise = 1000,ff = 10000,a = False):
    if fibers == 1:
        return opthext0(fwhm, npix, size, noise, ff, a)
    elif fibers == 7:
        return opthext1(fwhm, npix, size, noise, ff, a)
    elif fibers == 19:
        return opthext2(fwhm, npix, size, noise, ff, a)
    else:
        print("Please put in a valid number of fibers (1,7, or 19)")
        exit()

#make a rectangle of 1s in a square of 0s
def makerect(w, h):
    if w > h:
        rect = np.zeros((w,w))
        for i in range(w):
            if i >= (w-h)/2 and i < (w+h)/2:
                rect[i,:] = 1
    else:
        rect = np.zeros((h,h))
        for i in range(h):
            if i >= (h-w)/2 and i < (w+h)/2:
                rect[:,i] = 1
    return rect

#do optimal extraction on a slit with multiple pixels 
#pixelsize given in fraction of total slit length
def optrext(pixelsize, fwhm, npix=100, noise=50, ff=10000, a=False, m=True):
    #define variables
    width = .75 #width of slit in arcsec
    std = fwhm/2.355 #convert fwhm to stdev
    size = 1 #length of slit in arcsec
    pixheight = pixelsize*npix/size #height of pixel in terms of npix pixels

    #make a rectangular pixel mask to convolve over psf
    pixel = np.ones((int(pixelsize*npix/size), int(npix*width/size)))
    pixels = int(size/pixelsize)

    #make grid of positions across pixel space
    x = np.linspace(-.5*size,.5*size,npix)
    y = np.linspace(-.5*size,.5*size,npix)
    pos = [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]

    #define 2d pdf, evaluate to make a grid of values, normalize
    if m:
        tpsf = moffat((0,0), x, y, fwhm, 2.9)
        tpsf = tpsf/np.sum(tpsf)
    else:
        pdf = stats.multivariate_normal([0,0], [[std,0],[0,std]])
        tpsf = pdf.pdf(pos)

    #make sure each has appropriate sum
    tpsf = tpsf/np.sum(tpsf)
    psf = tpsf * ff

    #convolve fibers with psfs to get simulated fiber data
    #(sorry about how convoluted these lines are, I kind of was just having fun
    #with how much I could do at once by abusing lots of stuff)
    #'''
    pixellocs = np.asarray(list(zip(np.linspace(0,npix,pixels+1) + pixheight/2,
        np.zeros(pixels+1)+npix//2))).astype(int)
    data = addnoise(np.random.poisson(np.asarray([convolve(pixellocs[i], 
        pixel, psf) for i in range(pixels)])), noise)
    tdata = np.asarray([convolve(pixellocs[i], pixel, tpsf)
        for i in range(pixels)])
    '''
    alldata = convolve2d(psf, pixel, mode = 'valid')
    talldata = convolve2d(tpsf, pixel, mode = 'valid')
    ypos = list(range(0,len(alldata)-int(pixheight),int(pixheight)))\
            + [len(alldata)-1]
    data = addnoise(np.random.poisson(alldata[ypos]), noise)
    tdata = talldata[ypos]
    '''
    var = np.square(noise) + data #variance of noise
    return optimalextract(data, tdata, var, aperture = a) 

#return a coordinate grid with ordered pairs of each x,y combination
#equivalent to [[[x[i],y[j]] for j in range(npix)] for i in range(npix)]
#barely even faster, but there's probably a more elegant way to do this
def makegrid(x, y):
    points = np.asarray(np.meshgrid(x, y, indexing = 'ij'))
    return np.swapaxes(np.swapaxes(points, 0, 1), 1, 2)
