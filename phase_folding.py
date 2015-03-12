#
#   chuang
#   02.25.2015
#

"""
    Steps:
    
    Things to test: in addition to alpha, change the start date of the shift
"""

import numpy as np
from scipy import optimize
import scipy.stats
import sys
import matplotlib.pyplot as plt

import ft
import change_period as chper
import pdm
import period_fitting as pf

def ceph_fold(mjd, period, alpha, func=None, d0=0):
    """
        Takes in some dates and a starting period and then a period change of alpha. Folds the cepheids in just a little complicated manner
        
        Returns the phases.
    """
    phases = []
    #   zero the julian dates
    #mjd = mjd - min(mjd)
    if func == None:
        for i in mjd:
            period_new = period + alpha*(i+d0-min(mjd))
            phase = i/period_new - np.fix(i/period_new)
            phases.append(phase)
    else:
        for i in mjd:
            period_new = period + func((i+d0-min(mjd))/alpha)
            phase = i/period_new - np.fix(i/period_new)
            phases.append(phase)
    return phases

def fold_error_func(alpha, d0, mjd, mag_v, coeff, period, errors, order):
    new_phases = ceph_fold(mjd, period, alpha, d0=d0)
    template_magnitudes = []
    for i in new_phases:
        template_magnitudes.append(ft.cosine_series_phase(coeff, i, order))
    new_phases = np.asarray(new_phases)
    template_magnitudes = np.asarray(template_magnitudes)
    return np.asarray((template_magnitudes- mag_v)/errors)

def find_sec_change(mjd, mag_v, coeff, period, mag_v_err, order, num_alpha=201, alpha=0.0, d0=0):
    """
        Tests out a grid of values around an initial guess for alpha.
        
        Basically takes the mjd and the magnitudes of the cepheid, folds the data with various functions, and then for each resulting phase, compares it with the phase point using just the template.
    """
    alpha_range = np.linspace(-0.000005+alpha, 0.000005+alpha, num_alpha)
    #total_error = []
    num_points = len(mjd)
    red_chisq_n = num_points - 2
    #d0_range = np.linspace(-1000, 1000, 51)
    d0_range = [0]
    red_chisq = []
    for d in d0_range:
        print "Now trying d0 =", d
        for i in alpha_range:
            #print "Now trying alpha =", i
            i_errors = fold_error_func(i, d, mjd, mag_v, coeff, period, mag_v_err, order)
            #i_errors = np.asarray([b**2 for b in i_errors])
            #tot_error = np.sum(i_errors)
            tot_error = (np.array(i_errors)**2).sum()
            #total_error.append(tot_error)
            red_chisq.append(tot_error/red_chisq_n)
#   min_ind = np.argmin(total_error)
    min_ind = np.argmin(red_chisq)
    best_fit_alpha = alpha_range[min_ind]
    #best_d0 = d0_range[min_ind]
    return best_fit_alpha, red_chisq

def sec_change_per(mjd, mag_v, coeff, period, mag_v_err, order, num_periods=51,num_alpha=51, alpha=0.0, d0=0, a_range=0.000001, p_range=0.25):
    """
        Tests out a grid of values around an initial guess for alpha and period.
        
        Basically takes the mjd and the magnitudes of the cepheid, folds the data with various functions, and then for each resulting phase, compares it with the phase point using just the template.
    """
    alpha_range = np.linspace(-a_range+alpha, a_range+alpha, num_alpha)
    #total_error = []
    num_points = len(mjd)
    num_params = 2
    red_chisq_n = num_points - num_params # check this later
    per_range = np.linspace(-p_range + period, p_range + period, num_periods)
    # make an array of the appropriate size and populate it?
    red_chisq_arr = np.zeros((num_periods, num_alpha))
    for p in np.arange(len(per_range)):
        print "Now trying period =", per_range[p]
        for i in np.arange(len(alpha_range)):
            i_errors = fold_error_func(alpha_range[i], 0, mjd, mag_v, coeff, per_range[p], mag_v_err, order)
            tot_error = (np.array(i_errors)**2).sum()
            red_chisq_arr[p, i] = tot_error/red_chisq_n
    flat_ind = red_chisq_arr.argmin()
    bests = np.unravel_index(flat_ind, red_chisq_arr.shape)
    best_fit_period = per_range[bests[0]]
    best_fit_alpha = alpha_range[bests[1]]
    return best_fit_alpha, best_fit_period, red_chisq_arr

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

def best_chi(chi_square_array):
    """
        Stupid function for finding the best chi_square in an array of chi_squares
    """
    flat_midx = chi_square_array.argmin()
    best_idx = np.unravel_index(flat_midx, chi_square_array.shape)
    best_chi_square = chi_square_array[best_idx]
    return best_chi_square

def min_bin(data, nbins=20, min_elements=15, bin_min=0.0, bin_max=1.0, bins=None):
    """
       Returns the bin widths that bin the data, but combines bins that have fewer than 15 points (default) per bin.
    """
    if bins==None:
        bins = np.linspace(bin_min, bin_max, nbins+1)
        #print bins
    hist = np.digitize(data, bins)
    while len(np.where(np.bincount(hist)< min_elements)[0]) > 1:
        combine_bins = np.where(np.bincount(hist) < min_elements)[0]
        comb_bins = combine_bins[1:len(combine_bins)]
        for i in comb_bins:
            #print "number of bins", len(bins)
            #print "bin to delete", i
            bins = np.delete(bins, i)
            hist = np.digitize(data, bins)
            if np.where(np.bincount(hist)< min_elements) > 1:
                break
    return bins

def simple_sigma_clip(mag_v, phases, sigma=3, nbins=20, min_elements=15):
    """
        Bins the data, then clips points that fall more than 3 sigma off of the median value, as determined by the variance of the bin (default). Sets a minimum number of points per bin to be 15.
        
        Returns the indicies of the good data.
    """
    good_bins = min_bin(phases, nbins=nbins, bin_min=0.0, bin_max=1.0)
    binned_dat_medians = scipy.stats.binned_statistic(phases, mag_v, bins=good_bins, statistic='median')
    binned_dat_errors = scipy.stats.binned_statistic(phases, mag_v, bins=good_bins, statistic=np.std)
    medians = binned_dat_medians[0]
    sigmas = binned_dat_errors[0]
    cut = sigmas*sigma
    bin_numbers = binned_dat_medians[2]
    unique_bin_numbers = set(bin_numbers)
    good_inds = []
    for i in unique_bin_numbers:
        inds = np.where(bin_numbers == i)[0]
        median_pt = medians[i-1]
        sig_range = cut[i-1]
        for n in inds:
            if ((median_pt - sig_range) < mag_v[n]) & (mag_v[n] < (median_pt + sig_range)):
                good_inds.append(n)
    return np.asarray(good_inds)

def better_sigma_clip():
    print "TODO"



#####

#   plot everything doubled and see how well it agrees

if __name__=="__main__":
    cepheid = 'xycar'
    #cepheid = 'svul'
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

    #   make a template cepheid
    coeff, order = ft.fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=7)
    #per, theta = pdm.find_period_arb(new_mjds, good_mag_v, period)
    #per, theta = pf.PDM_precise(new_mjds, good_mag_v, period)
    
    #   create cepheid data with some complex phase.
    shifted_mags = ft.evolving_cosine(coeff, good_mjd, 4, period, alpha, ft.linear_per_change, min_date = min(mjd))
    shifted_phases = good_mjd/period - np.fix(good_mjd/period)
    test_phases = np.asarray(ceph_fold(good_mjd, period, alpha))
    
    #   add in some Gaussian noise
    noise = np.random.normal(0, 0.03, len(shifted_mags))
    mags_with_noise = shifted_mags + noise

    #   can't really make an o-c diagram now because instead of shifting the dates as I was before, I am now just shifting the magnitudes.
    
    phase_points = np.linspace(0, 1.0, 201)
    template_mags = ft.cosine_series_phase(coeff, phase_points, order)

    #   plot to make sure that it's actually doing that
    plt.scatter(shifted_phases, mags_with_noise, color='red', alpha=0.3)
    #plt.scatter(good_phases, good_mag_v, color='red', alpha=0.3)
    #plt.scatter(good_phases, good_mag_v, color='grey', alpha=0.3)
    plt.scatter(test_phases, shifted_mags, color='blue', alpha=0.3)
    #plt.scatter(phase_points, template_mags, color='blue', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.xlabel("Phase")
    plt.ylabel("V-Band Magnitude")
    plt.savefig("../output/{}_real_data.png".format(cepheid))
    plt.show()
    
    
    num_alpha = 201
    #sys.exit()
    #   find the change that works the best
    #best_alpha, chi_squares = find_sec_change(good_mjd, mags_with_noise, coeff, period, good_mag_v_err, order)
    best_alpha, chi_squares = find_sec_change(good_mjd, good_mag_v, coeff, period, good_mag_v_err, order, num_alpha=num_alpha)
    best_phases = np.asarray(ceph_fold(good_mjd, period, best_alpha))
    temp_mags = ft.cosine_series_phase(coeff, best_phases, order)
    temp_phases = np.asarray(ceph_fold(good_mjd, period, best_alpha))
    alpha_range = np.linspace(-0.000005+alpha, 0.000005+alpha, num_alpha)
    
    plt.scatter(alpha_range, chi_squares, color='grey', alpha=0.3)
    plt.xlim([-0.000006, 0.000006])
    plt.ylabel("Chi-Square")
    plt.xlabel("Alpha")
    plt.savefig("../output/{}_real_chi_squares.png".format(cepheid))
    plt.show()
    
    #sys.exit()

    #   repeat once for the plotting
    x, y = duplicate(best_phases, mags_with_noise)
    #x, y = duplicate(best_phases, good_mag_v)
    #xor, yor = duplicate(temp_phases, temp_mags)
    #xor, yor = duplicate(good_phases, good_mag_v)
    xor, yor = duplicate(shifted_phases, mags_with_noise)

    #plt.scatter(x, y, color='grey', alpha=0.3)
    plt.scatter(best_phases, good_mag_v, color='grey', alpha=0.3)
    plt.scatter(good_phases, good_mag_v, color='red', alpha=0.3)
    #plt.scatter(xor, yor, color='blue', alpha=0.3)
    #plt.scatter(xor, yor, color='red', alpha=0.3)
    #plt.scatter(test_phases, shifted_mags, color='blue', alpha=0.3)
    plt.gca().invert_yaxis()
    plt.ylabel("V-Band Magnitudes")
    plt.xlabel("Phase")
    plt.savefig("../output/{}_fake_final.png".format(cepheid))
    plt.show()



