#
#   chuang
#   02.25.2015
#

"""
    Steps:
    
    1. Use pdm analysis to find the best fit period.--at the moment the fourier analysis problem just uses the literature period though
    2. May want to do a three-sigma clip on the real data, if you are using real data. 
    3. Create a template using the best fit period and the Cepheid data. 
    4. Take real data and fold with a changing period.
    5. With the output phases from this folding, compare with the magnitude that the Cepheid would be at for that part of its phase (based off of the template).
    6. Minimize with weighted least squares, but trying a grid of different alphas, etc.
    7. Test and see if you can recover your result.
"""

import numpy as np
from scipy import optimize
import sys
import matplotlib.pyplot as plt

import ft
import change_period as chper
import pdm
import period_fitting as pf

def ceph_fold(mjd, period, alpha):
    """
        Takes in some dates and a starting period and then a period change of alpha. Folds the cepheids in just a little complicated manner
        
        Returns the phases.
        """
    phases = []
    #   zero the julian dates
    #mjd = mjd - min(mjd)
    for i in mjd:
        period_new = period + alpha*(i-min(mjd))
        phase = i/period_new - np.fix(i/period_new)
        phases.append(phase)
    return phases

def fold_error_func(alpha, mjd, mag_v, coeff, period, errors, order):
    new_phases = ceph_fold(mjd, period, alpha)
    template_magnitudes = []
    for i in new_phases:
        template_magnitudes.append(ft.cosine_series_phase(coeff, i, order))
    new_phases = np.asarray(new_phases)
    template_magnitudes = np.asarray(template_magnitudes)
    return np.asarray((template_magnitudes- mag_v)/errors)

def find_sec_change(mjd, mag_v, coeff, period, mag_v_err, order, alpha=0.0):
    """
        Tests out a grid of values around an initial guess for alpha.
        
        Basically takes the mjd and the magnitudes of the cepheid, folds the data with various functions, and then for each resulting phase, compares it with the phase point using just the template.
    """
    alpha_range = np.linspace(-0.000005+alpha, 0.000005+alpha, 201)
    total_error = []
    for i in alpha_range:
        print "Now trying alpha =", i
        i_errors = fold_error_func(i, mjd, mag_v, coeff, period, mag_v_err, order)
        i_errors = np.asarray([b**2 for b in i_errors])
        tot_error = np.sum(i_errors)
        total_error.append(tot_error)
    min_ind = np.argmin(total_error)
    best_fit_alpha = alpha_range[min_ind]
    return best_fit_alpha

def find_sec_chng(phases, good_mjd, mag_v, coeff, period, mag_v_err, order, alpha=0.0):
    """
        It would be better if this were given in terms of phases...
        
        So the phases are going to be the phases of the changing Cepheid (so they will just be its mjd/period, where period is constant). Then the template is going to be divided by changing periods until they two look sufficiently similar.
        
        NOT DONE
    """
    alpha_range = np.linspace(-0.000005+alpha, 0.000005+alpha, 101)
    total_error = []
    for i in alpha_range:
        print "trying alpha =", i
        i_errors = fold_error_func(i, mjd, mag_v, coeff, period, mag_v_err, order)
        i_errors = np.asarray([b**2 for b in i_errors])
        tot_error = np.sum(i_errors)
        total_error.append(tot_error)
    min_ind = np.argmin(total_error)
    best_fit_alpha = alpha_range[min_ind]
    return best_fit_alpha

def duplicate(xdata, ydata, num=2):
    i = 1
    total_xdata = np.asarray(xdata)
    total_ydata = np.asarray(ydata)
    while i < 2:
        new_data = xdata + 1.0*i
        total_xdata = np.concatenate([total_xdata, new_data])
        total_ydata = np.concatenate([total_ydata, ydata])
        i = i + 1
    return total_xdata, total_ydata


def unfold_simple(phases, mjd, start_period, alpha):
    """
        Basically does the opposite of ceph_fold. It takes Cepheid data that's been phased by a changing period and then it 'unfolds' it by the simple linear model.
        
        I realized that all this has to do is multiply the Cepheid phases by the fixed period..., but hopefully this will be a bit cleaner since you can account for the zeropointing.
        
        This doesn't work since the data is periodic. I mean the only way it could is if I included the epoch data, but is that a good idea?
    """
    n_el = len(phases)
    epoch = np.fix(mjd/start_period)
    mjd = []
    for i in np.arange(n_el):
        mjd_el = (phases[i]*start_period)/(1.-phases[i]*alpha) + epoch[i]*start_period
        mjd.append(mjd_el)
    mjd = np.asarray(mjd)
    return mjd


#####

#   Suggestion from Kevin to make fake data--just change the frequency in the cosine function to have it depend on time like the function you want. Then you have template "fake" data.

#   plot everything doubled and see how well it agrees

if __name__=="__main__":
    cepheid = 'xycar'
    ceph_data = ft.get_sorted_data(cepheid)
    mjd = ceph_data[:,0]
    mag_v = ceph_data[:,1]
    epochs = ceph_data[:,2]
    phases = ceph_data[:,3]
    mag_v_err = ceph_data[:,4]

    #   do a crappy clipping of the really bad points--fix this later
    good_ind = np.where(mag_v < 10)
    good_ceph_dat = ceph_data[good_ind]
    good_mjd = mjd[good_ind]
    good_mag_v = mag_v[good_ind]
    good_epochs = epochs[good_ind]
    good_phases = phases[good_ind]
    good_mag_v_err = mag_v_err[good_ind]
    
    period = ft.get_period(cepheid)
    alpha = 0.0000001
    
    #   change the period of the cepheid
    #new_mjds, new_periods = chper.secular_change(good_mjd, period, 0.001)

    #   make a template cepheid
    #coeff, order = ft.fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=7)
    #per, theta = pdm.find_period_arb(new_mjds, good_mag_v, period)
    #per, theta = pf.PDM_precise(new_mjds, good_mag_v, period)
    
    #   create template "real" cepheid data
    coeff, order = ft.fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=7)
    shifted_mags = ft.evolving_cosine(coeff, good_mjd, 7, period, alpha, ft.linear_per_change, min_date = min(mjd))
    
    shifted_phases = good_mjd/period - np.fix(good_mjd/period)
    test_phases = np.asarray(ceph_fold(good_mjd, period, alpha))

    #new_phases = ceph_fold(good_mjd, period, alpha)
    #new_phases = new_mjds/per - np.fix(new_mjds/per)
    #new_mjds = unfold_simple(new_phases, mjd, period, alpha)
    
    #temp_mags = ft.cosine_series(coeff, good_mjd, order, period)
    #temp_st_phs = good_mjd/period - np.fix(good_mjd/period)
    
    #   plot to make sure that it's actually doing that
    plt.scatter(shifted_phases, shifted_mags, color='red', alpha=0.3)
    plt.scatter(good_phases, good_mag_v, color='grey', alpha=0.3)
    plt.scatter(test_phases, shifted_mags, color='blue', alpha=0.3)
    #plt.gca().invert_yaxis()
    plt.show()
    
    #sys.exit()
    #test_phases = np.linspace(0, 1.0, 1000)
    #template_mags = ft.cosine_series_phase(coeff, test_phases, order)
    
    #plt.scatter(new_phases, good_mag_v, color='red', alpha=0.3)
    #plt.scatter(test_phases, template_mags, color='grey', alpha=0.3)
    #plt.gca().invert_yaxis()
    #plt.show()
    
    #   get the period
    #period = ft.get_period(cepheid)

    #   find the change that works the best
    best_alpha = find_sec_change(good_mjd, shifted_mags, coeff, period, good_mag_v_err, order)
    #best_alpha = 0.0

    best_phases = np.asarray(ceph_fold(good_mjd, period, best_alpha))
    

    #   template magnitudes
    template_mags = ft.cosine_series_phase(coeff, best_phases, order)

    x, y = duplicate(best_phases, shifted_mags)
    #xor, yor = duplicate(best_phases, template_mags)

    #phases_real = good_mjd/period - np.fix(good_mjd/period)

    plt.scatter(x, y, color='red', alpha=0.3)
    plt.scatter(xor, yor, color='grey', alpha=0.3)
    plt.scatter(test_phases, shifted_mags, color='blue', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.show()



