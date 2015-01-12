#	chuang
#	12.10.14
#
#	module that fits some cepheid data with the appropriate periods. Options to use various methods based on what the data looks like (for now we'll just consider PDM=phase dispersion minimization from PyTiming)
#

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from PyAstronomy.pyTiming import pyPDM
import cepheid_plotting
from scipy import interpolate
from scipy.interpolate import UnivariateSpline

#   Fits period to data using phase dispersion minimization, and plots the theta values. Good for data that is sparse and has a lot baseline. This is also good when you have a decent initial guess to the period of the Cepheid (within a 4-day range)
#   Need to get rid of the odd phase shift you have.
#
#   Parameters:
#   t = the dates
#   y = the magnitudes
#   period = best guess for the period
#   precision = how precise you want the period to be, 0.001 if nothing given

def extract_median(t, y):
    binned_t = []
    binned_mag = []
    t_diff = np.diff(t)
    end_ind = np.concatenate([np.where(t_diff > 0.04)[0] + 1, [len(t)]], 1)
    begin_ind = np.concatenate([[0], np.where(t_diff > 0.04)[0] + 1],1)
    bin_number = np.arange(0, len(end_ind))
    for i in bin_number:
        bin = t[begin_ind[i]:end_ind[i]]
        bin_mag = y[begin_ind[i]:end_ind[i]]
        t_median = np.median(bin)
        m_median = np.median(bin_mag)
        binned_t.append(t_median)
        binned_mag.append(m_median)

    binned_t = np.asarray(binned_t)
    binned_mag = np.asarray(binned_mag)
    return binned_t, binned_mag


def PDM_precise(t, y, period, name='none', precision=0.001):

    minperiod = period - 2.0
    maxperiod = period + 2.0
    S = pyPDM.Scanner(minVal=minperiod, maxVal=maxperiod, dVal=precision, mode="period")
    P = pyPDM.PyPDM(t, y)
    print "Trying fit..."
    f1, t1 = P.pdmEquiBinCover(10, 3, S)
    f2, t2 = P.pdmEquiBin(10, S)

    best_fit_idx = t1.argmin()
    revised_period = f1[best_fit_idx]
    print "Best fit period:", revised_period
    
    phase_shift = y[0]/revised_period - y[0]/period
    
    pp, yy = cepheid_plotting.plot_array(t, y, revised_period)
    pp2, yy2 = cepheid_plotting.plot_array(t, y, period)

#    pp_ind = np.argsort(pp)
#    pp = pp[pp_ind]
#    yy = yy[pp_ind]
    
#    deltaphases = np.linspace(0, 1, 1001)
#    s = UnivariateSpline(pp2, yy2)
#    phasemags = s(deltaphases)
#    phaseshift = deltaphases[np.argmax(phasemags)]
#    print phaseshift
#    pp = pp - phaseshift

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.xlabel("Frequency")
    plt.ylabel("Theta")
    plt.plot(f1, t1, 'bp-')
    plt.plot(f2, t2, 'gp-')
    plt.legend(["pdmEquiBinCover", "pdmEquiBin"])
    plt.text(minperiod, 1.1, "best fit period = {} days".format(revised_period), fontsize=14)
    plt.subplot(2, 1, 2)
    plt.scatter(pp, yy, alpha=0.7, color='red')
    #plt.scatter(pp2, yy2, alpha=0.7, color='grey')
    plt.gca().invert_yaxis()
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Phase")
    if name!='none':
        plt.savefig("../output/{}_pdm_fitting.png".format(name))
        plt.title("Result of PDM analysis for {}".format(name))
    else:
        plt.title("Result of PDM analysis")
    plt.show()

#   Input parameter depends on if the period is changing or not.
#   If the period is not changing, then use all of the data or the binned data to fit
#   If the period IS changing, then fit the spline to a specific place.
#
#   Parameters:
#   p = phase to get a value in days (just one period)
#   y = the magnitudes

def lc_shape(p, y, period):

    d = p*period
    s = interpolate.InterpolatedUnivariateSpline(d, y)
    return s

