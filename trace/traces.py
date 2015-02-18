# coding: utf-8

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import math

from astropy.io import fits

from peakdetection import peak_detection_mean_window

def wcs_to_pix(x):
    return int(math.floor(x + 0.5))

def linear_fit(xx, yy):
    xx = np.asarray(xx)
    yy = np.asarray(yy)
    xm = xx.mean()
    ym = yy.mean()
    xn = xx - xm
    yn = yy - ym
    up = (xn * yn).sum()
    down = (xn * xn).sum()
    bet = up / down
    alp = ym - bet * xm
    return bet, alp

def prediction(xx, yy, x, n=5):

    # take n points
    nn = min(len(xx), len(yy))
    if nn == 0:
        return 0.0
    elif nn == 1:
        return yy[0]
    else:
        # fit with nn points
        xf = xx[:nn]
        yf = yy[:nn]
        bet, alp = linear_fit(xf, yf)
        return bet * x + alp
    
    return 0.0

def fit_para_equal_spaced(dd):
    '''Parabole that passes through 3 points
        
    With X=[-1,0,1]
    '''
    
    Y2 = dd[1]
    Y1 = dd[0]
    Y3 = dd[2]

    C = Y2        
    B = 0.5 * (Y3 - Y1)
    A = 0.5 * (Y1 + Y3 - 2 * Y2)
    return A, B, C


def interp_max_3(dd):
    '''Parabola that passes through 3 points
        
    With X=[-1,0,1]
    '''

    A, B, C = fit_para_equal_spaced(dd)
    
    return -B / (2*A), C - B * B / (4*A)


def delicate_centre(x, y):
    pos = np.polyfit(x, y, deg=2)
    tx = -pos[1] / (2 * pos[0])
    py = np.polyval(pos, tx)
    return tx, py, pos

TRACE_BACKWARD, TRACE_FORWARD = 0, 1

class Trace(object):
    def __init__(self, start=0, direction=TRACE_FORWARD, npred=5):
        self.start = start
        self.lost = None
        self.sample_c = []
        self.trace_c = []
        self.peak_c = []
        self.sample_f = []
        self.trace_f = []
        self.peak_f = []
        self.npred = npred
        self.direction = direction

        self.set_direction(direction)
        
    def set_direction(self, direction):
        if direction != self.direction:
            self.reverse()
            self.direction = direction
            if self.direction == TRACE_BACKWARD:
                self.predict = self.predict_backward
            elif self.direction == TRACE_FORWARD:
                    self.predict = self.predict_forward
            else:
                raise ValueError('direction is not F or B')

    def predict_forward(self, col):

        xpred = self.sample_f[-self.npred:]
        ypred = self.trace_f[-self.npred:]

        expected = prediction(xpred, ypred, col, n=self.npred)
        return expected

    def predict_backward(self, col):
        xpred = self.sample_f[:-(self.npred+1):-1]
        ypred = self.trace_f[:-(self.npred+1):-1]

        expected = prediction(xpred, ypred, col, n=self.npred)
        return expected

    def reverse(self):
        self.sample_c.reverse()
        self.trace_c.reverse()
        self.peak_c.reverse()
        self.sample_f.reverse()
        self.trace_f.reverse()
        self.peak_f.reverse()

#
# Multiple inheritance here?
#
class FiberTrace(Trace):
    def __init__(self, fibid, boxid, start=0):
        super(FiberTrace, self).__init__(start=start)
        self.boxid = boxid
        self.fibid = fibid

    def __str__(self):
        return "FiberTrace(fibid=%i, boxid=%i, start=%i)" % (self.fibid, self.boxid, self.start)


def init_traces(image, center, hs, background, npred, maxdis=9.0):

    ixmin= 100
    ixmax = 4000

    xx = np.arange(image.shape[0])


    cut_region = slice(center-hs, center+hs)
    cut = image[:,cut_region]
    colcut = cut.mean(axis=1)
    maxt = peak_detection_mean_window(colcut, x=xx, k=3, xmin=ixmin, xmax=ixmax, background=background)
    npeaks = len(maxt)
    peakdist = np.diff(maxt[:,1])
    # number of peaks
    print('npeaks', npeaks)
    print('distances between peaks (max)', peakdist.max(), '(min)', peakdist.min())
    #
    fiber_traces = {}
    fibid = 1
    gcounter = 0
    boxid = 1
    for xpeak, dis in zip(maxt[:,1], peakdist):
        gcounter += 1
        if dis > maxdis:
            print(xpeak, dis, gcounter, boxid, fibid)
            gcounter = 0
            boxid += 1
        fiber_traces[fibid] = FiberTrace(fibid, boxid, start=center)
        fibid += 1

    fw = 2

    print("initialice the centers of fibers")
    for fibid, trace in fiber_traces.items():
        trace.sample_c.append(center)
        trace.trace_c.append(maxt[fibid,1])
        trace.peak_c.append(maxt[fibid,2])
        pixmax = int(maxt[fibid,0])
        # Take 2*2+1 pix
        # This part and interp_max_3(image[nearp3-1:nearp3+2, col])
        # should do the same
        tx, py, pos = delicate_centre(xx[pixmax-fw: pixmax+fw+1], 
                                     colcut[pixmax-fw: pixmax+fw+1])
        trace.sample_f.append(center)       
        trace.trace_f.append(tx)
        trace.peak_f.append(py)
    print('Done')
    return fiber_traces

def we_are_in_limits(col):
    return (col > 6) and (col < 4090)

print('Loop over traces')

maxdis = 2.0

def trace_backward(trace, image, start, hs, background):
    return trace_common(trace, image, start, hs, direction=TRACE_BACKWARD, background=background)
    
def trace_forward(trace, image, start, hs, background):
    return trace_common(trace, image, start, hs, direction=TRACE_FORWARD, background=background)
    
def trace_common(trace, image, start, hs, direction, background):

    if direction == TRACE_FORWARD:
        dirmod = 1
    elif direction == TRACE_BACKWARD:
        dirmod = -1
    else:
        raise ValueError('direction must be either F or B')

    trace.set_direction(direction)

    col = start
    while we_are_in_limits(col):
        print('we are in column', col)    
        col = col + dirmod * hs
        print('we go to column', col)

        expected = trace.predict(col)
        print('we predict the peak will be in coordinate', expected)

        epix = wcs_to_pix(expected)

        print('extract a region around the expected peak')
        region = image[epix-5:epix+5,col-hs:col+hs].mean(axis=1)
        thisx = np.arange(epix-5,epix+5)
        print('find the peak')
        maxt = peak_detection_mean_window(region, x=thisx, background=background)
        print(maxt)
        
        if len(maxt) < 1:
            print('no peaks, exit, col is ', col)
#            thisx = np.arange(epix-15,epix+15)
#            region = image[epix-15:epix+15,col-hs:col+hs].mean(axis=1)
#            plt.plot(thisx, region, 'r*-')
        
#            plt.show()
            break


        nearp1 = np.argmin(np.abs(maxt[:,1] - expected))
        nearp2 = maxt[nearp1,1]
        nearp3 = int(nearp2)

        print('check the peak is not further than npixels')
        if abs(nearp2 - expected) > maxdis:
            print('peak too far from prediction')
            break
        else:
            print('correct')

        print('fit the peak', nearp3)
        xthin, ythin = interp_max_3(image[nearp3-1:nearp3+2, col])
        xthin += nearp3
        print('fit the peak', xthin, ythin)

        #plt.xlim([235,245])
        #plt.plot(thisx, region, 'r*-')
        #plt.plot(maxt[:,1], 1.1 * maxt[:,2], 'g*')
        #plt.plot([xthin], [1.1 * ythin], 'bo')
        #plt.show()


        trace.sample_f.append(col)
        trace.trace_f.append(xthin)
        trace.peak_f.append(ythin)

        trace.sample_c.append(col)
        trace.trace_c.append(nearp2)
        trace.peak_c.append(maxt[nearp1,2])

        print('go to the next point')
    return trace


if __name__ == '__main__':
    print("load image")
    cmap = plt.cm.get_cmap('afmhot')
    image1 = fits.getdata('product_002.fits')

    print('trace detection')

    cstart = 2000
    hs = 5
    background1 = 150.0
    npred = 5

    fiber_traces = init_traces(image1, center=cstart, hs=hs, 
                                background=background1, npred=npred)

    i = 1
    for trace in fiber_traces.values():

        col = trace.sample_c[0]
    
        trace = trace_backward(trace, image1, cstart, hs, background1)

        trace = trace_forward(trace, image1, cstart, hs, background1)

        plt.plot(trace.sample_f, trace.trace_f, 'r*')
        plt.savefig('fig%i.png' % i)
        i += 1
    
    # check if samples are sorted
    #assert all(a <= b for a, b in zip(trace.sample_f, trace.sample_f[1:]))
    

    

