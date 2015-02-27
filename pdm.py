#
#   chuang
#   02.20.15
#
#

import numpy as np
import matplotlib.pylab as mpl
from PyAstronomy.pyTiming import pyPDM
import asciitable

"""
    Code to perform phase-dispersion analysis to find an optimal period.
"""

def find_period(lc_file, guess_period, stepsize=0.01, bins=15):
    """
        Finds the period of a Cepheid given an initial guess (use the literature data for this) and then tries to find a better period by minimizing the dispersion. 
        
        Theta is the goodness of fit. Set some tolerance for theta, so that we have a cutoff about what is changing period and what is not changing period.
    """
    lc_dat = asciitable.read(lc_file)
    y = lc_dat['V-band magnitude']
    x = lc_dat['mjd']
    S = pyPDM.Scanner(minVal=guess_period-0.75, maxVal=guess_period+0.75, dVal=stepsize, mode="period")
    P = pyPDM.PyPDM(x, y)
    f1, t1 = P.pdmEquiBinCover(15, 3, S)
    min_theta = min(t1)
    min_ind = np.argmin(min_theta)
    period = f1[min_ind]
    return period, min_theta

def find_period_arb(mjd, mag, guess_period, stepsize=0.01, bins=15):
    """
        Finds the period of a Cepheid given an initial guess (use the literature data for this) and then tries to find a better period by minimizing the dispersion.
        
        Theta is the goodness of fit. Set some tolerance for theta, so that we have a cutoff about what is changing period and what is not changing period.
    """
    y = mag
    x = mjd
    S = pyPDM.Scanner(minVal=guess_period-0.5, maxVal=guess_period+0.5, dVal=stepsize, mode="period")
    P = pyPDM.PyPDM(x, y)
    f1, t1 = P.pdmEquiBinCover(15, 3, S)
    min_theta = min(t1)
    min_ind = np.argmin(min_theta)
    period = f1[min_ind]
    return period, min_theta
