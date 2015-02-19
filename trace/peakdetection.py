import numpy as np
from scipy.ndimage.filters import generic_filter

def _vecS1(k, data):
    '''max filter for peak detection.'''

    def func(x):
        return x[k] - 0.5 * (x[:k].max() + x[k+1:].max())
    
    ap = generic_filter(data, func, size=2*k+1)
    return ap

def _vecS2(k, data):
    '''min filter for peak detection.'''
    
    def func(x):
        return x[k] - 0.5 * (x[:k].mean() + x[k+1:].mean())
    
    ap = generic_filter(data, func, size=2*k+1)
    return ap

# http://tcs-trddc.com/trddc_website/pdf/srl/palshikar_sapdts_2009.pdf

def _generic_peak_filtering(method, y, x=None, k=3, h=1.0, background=0.0,
                            xmin=None, xmax=None):
    '''Helper class for filtered peak finding.'''
    
    if x is None:
        x = np.arange(len(y))

    if xmin is None:
        xmin = x[0]

    if xmax is None:
        xmax = x[-1]
    
    candidates = []
    # Fixme...
    
    ap = method(k, y)
    # mean and std of pos values
    appos = ap[ap>0]
    #
    # This does not work, we need something local, not global
    #mpos = appos.mean()
    #spos = appos.std()
    
    for idx, a in enumerate(ap):
        #if a > background and (a - mpos) > (h * spos):
        if a > 0 and y[idx] > background and (xmin <= x[idx] <= xmax):
            candidates.append(idx)

    # Filter close elements
    result = []
    while len(candidates) > 1:
        
        x0 = candidates[0]
        for idx, xi in enumerate(candidates[1:], 1):
            if abs(xi - x0) > k:
                break

        group = candidates[:idx]

        remaining = candidates[idx:]
        maxvalid = y[group].argmax()
        finalid = group[maxvalid]
        result.append((finalid, x[finalid], y[finalid]))
        
        candidates = remaining
    
    return np.array(result)

def peak_detection_mean_window(y, x=None, k=3, h=1.0, background=0.0, xmin=None, xmax=None):
    '''Detec peaks using a mean filter in a 2*k +1 window.

    The filter value is the central point minus the mean of the 2*k
    adjacent points.
    '''
    return _generic_peak_filtering(_vecS2, y, x=x, k=k, h=h,
                                   background=background, xmin=xmin, xmax=xmax)

def peak_detection_max_window(y, x=None, k=3, h=1.0, background=0.0, xmin=None, xmax=None):
    '''Detec peaks using a max filter in a 2*k +1 window.

    The filter value is the central point minus the mean of the max
    of the k-left and the max of the k-rigth neighbours.
    '''
    return _generic_peak_filtering(_vecS1, y, x=x, k=k, h=h,
                                   background=background, xmin=xmin, xmax=xmax)

