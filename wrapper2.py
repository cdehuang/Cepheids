#   chuang
#   1.22.15
#
#   Revised version of the wrapper code. Had to make too many changes so I just decided to start a new one.
#   This code should run the fourier transforming stuff and then produce a plot of chi-square.

import numpy as np
import matplotlib.pyplot as plt
import ft
from scipy import stats
from PyAstronomy.pyTiming import pyPDM
import matplotlib.pylab as mpl

#   Calculate the residuals

def residuals(observed, expected):
    return observed - expected

#   Calculate chi-squared --> the errors need to be assessed in a different way, but for now just get the errors on those specific data points I suppose?

def chi_square(resid, errors):
    cs = np.sum(np.power(resid,2.)/np.power(errors,2.))
    return cs

def chi_square2(resid, expected):
    cs = np.sum(np.power(resid,2.)/np.power(expected, 2.))
    return cs

#   not super useful right now, but basically it looks at the fit. it is just graphs at the moment, but you can get it to print stuff later too.

def examine_fit(name):
    # name = 'xycar' # <-- this should be the sole input parameter for now. later on the other input parameter will be the date of observation
    cepheid_data = ft.get_ceph_data(name, mag_shift=0.0)
    mjd = cepheid_data[:,0]
    mag_v = cepheid_data[:,1]
    epoch = cepheid_data[:,2]
    phase = cepheid_data[:,3]
    err = cepheid_data[:,4]
    period = ft.get_period(name)
    #   fold the data
    folded_phase, mags, errors = ft.fold_data(mjd, mag_v, period, err=err)
    #   get the medians of the folded data
    ph_med, run_med, err_med = ft.get_median(folded_phase, mags, errors=errors)
    #   get the actual fourier fit parameters, note that the input phase parameters were all in days
    coeff, ord = ft.fourier_series(ph_med*period, run_med, name)

    #   The fourier fitting program really struggles with trying to fit all of the data points.
    #   Try to get an idea about the quality of the fit by finding the chi square. First take the fit that you have and see what the chi-square is when you fit it with all of the data.

    expected = ft.cosine_series(coeff, folded_phase*period, ord, period)
    observed = mags

    resid = residuals(observed, expected)

    chisquare = chi_square2(resid, expected)

    dfolded_phase, dmags, derrors = ft.fold_data(mjd, mag_v, period, err=err, num=4)

    plt.scatter(folded_phase*period, observed, color='grey', alpha=0.2)
    plt.scatter(folded_phase*period, expected, color='red', alpha=0.3)
    #   plt.scatter(dfolded_phase*period, ft.cosine_series(coeff, dfolded_phase*period, ord, period), color='red', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.xlabel("Days")
    plt.ylabel("V-band magnitudes")
    plt.title("{}".format(name))
    plt.show()

    #   try plotting and then adding in like a linear change etc. and when you plot the whole thing fold all the data points onto the plot again so that you can see if it matches.

    phase_shifts = np.linspace(-0.5, 0.5, 101)
    chi_sq_array = []
    for i in phase_shifts:
        shifted_phases = (folded_phase+i)*period

        observed = mags
        expected = ft.cosine_series(coeff, shifted_phases, ord, period)

        resid = residuals(observed, expected)

        output = chi_square2(resid, errors) # <--need to decide which is the appropriate thing to use
        chi_sq_array.append(output)

    """
    plt.scatter(folded_phase*period, observed, color = 'grey', alpha=0.2)
    plt.scatter(folded_phase*period+1, expected, color='red', alpha=0.3)
    plt.show()
    """
    plt.scatter(phase_shifts, chi_sq_array)
    plt.ylabel("$\chi^2$")
    plt.xlabel("Phase Shift")
    plt.title("{}".format(name))
    plt.show()

#   constant period change <--adds a period change to some data, where the period is changing by a constant amount per year...i would do per cycle, but that seems like it would be hard to make consistent since the periods differ by so much between cepheid. Change is given in days
#   start off with a simple example, where you assume you start off with the base period and are just adding more periods onto it.
#   Assuming 365.25 days = 1 year

def constant_change(mjd, change, d0=0):
    change_per_day = change/365.25
    if d0 == 0:
        d0 = min(mjd)
    diff = np.asarray([mjd[i] - d0 for i in range(len(mjd))])
    adjusted_mjd = mjd + diff*change_per_day
    return adjusted_mjd

def cc_function(change, d0, mjd)
    change_per_day = change/365.25
    return

#   Try making something with alternative ways of period changes.

def nc_change():
    print "work in progress"

#   pdm analysis

def pdm_period(mjd, mags, period, dVal=0.01):
    minperiod = period - 1.0
    maxperiod = period + 1.0
    S = pyPDM.Scanner(minVal=minperiod, maxVal=maxperiod, dVal=dVal, mode="period")
    P = pyPDM.PyPDM(mjd, mags)
    f1, t1 = P.pdmEquiBinCover(20, 3, S)
    f2, t2 = P.pdmEquiBin(20, S)
    
    min_ind = np.argmin(t1)
    
    mpl.figure(facecolor='white')
    mpl.title("Result of PDM analysis")
    mpl.xlabel("Frequency")
    mpl.ylabel("Theta")
    mpl.plot(f1, t1, 'bp-')
    mpl.plot(f2, t2, 'gp-')
    mpl.legend(["pdmEquiBinCover", "pdmEquiBin"])
    mpl.show()
    return f1[min_ind]

#   try making some plots. Try making some plots where the period is changing very very slowly for xycar and then plot it against the actual data (after it's been folded)

name = 'xycar' # <-- this should be the sole input parameter for now. later on the other input parameter will be the date of observation
cepheid_data = ft.get_ceph_data(name, mag_shift=0.0)
mjd = cepheid_data[:,0]
mag_v = cepheid_data[:,1]
epoch = cepheid_data[:,2]
phase = cepheid_data[:,3]
err = cepheid_data[:,4]
period = ft.get_period(name)

adj_mjd = constant_change(mjd, 1)

folded_adj_phase, mags, errors = ft.fold_data(adj_mjd, mag_v, period, err=err)
folded_phase, mags, errors = ft.fold_data(mjd, mag_v, period, err=err)

plt.scatter(folded_adj_phase*period, mags, color='red', alpha=0.3)
plt.scatter(folded_phase*period, mags, color='grey', alpha=0.3)
plt.gca().invert_yaxis()
plt.title("{}".format(name))
plt.ylabel("V-band Magnitude")
plt.xlabel("Days")
plt.show()

#   Now run this through the dispersion-minimizing program and see what you get.

best_fit_period = pdm_period(adj_mjd, mags, period)
print best_fit_period

#   Now try to plot everything again with this new period

refit_phases, mags, errors = ft.fold_data(adj_mjd, mag_v, best_fit_period, err=err)
plt.scatter(folded_phase, mags, color='grey', alpha=0.3)
plt.scatter(folded_adj_phase, mags, color='red', alpha=0.3)
plt.scatter(refit_phases, mags, color='blue', alpha=0.3)
plt.gca().invert_yaxis()
plt.show()

#   If you think something has a constant change you can adjust the parameters by playing around with the 'd0' and with the change per year to minimize it. Try it with the actual xycar parameters


