#
#	chuang
#	03.09.2015
#

"""
    A script that makes a Monte Carlo given some fake data distributions.
    
    Compare the chi square, the output fit parameters, and the input parameters and see how well they match up.
"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii

import ft
import phase_folding as pf

#   number of iterations
n = 100

#   cepheid to base the data off of
cepheid = 'xycar'
real_period = ft.get_period('xycar')

ceph_data = ft.get_sorted_data(cepheid)
mjd = ceph_data[:,0]
mag_v = ceph_data[:,1]
epochs = ceph_data[:,2]
phases = ceph_data[:,3]
mag_v_err = ceph_data[:,4]

good_ind = pf.simple_sigma_clip(mag_v, phases)
good_ceph_dat = ceph_data[good_ind]
good_mjd = mjd[good_ind]
good_mag_v = mag_v[good_ind]
good_epochs = epochs[good_ind]
good_phases = phases[good_ind]
good_mag_v_err = mag_v_err[good_ind]

coeff, order = ft.fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=7)

#   starting parameter ranges, give these period and alpha uniformly distributed in those ranges
per_range = 0.25 # +/- from 'real' period, in days
alpha_range = 0.0000001 # +/-, also in days, this is going to be centered around 0
periods = np.random.uniform(real_period - per_range, real_period + per_range, n)
alphas = np.random.uniform(-alpha_range, alpha_range, n)
start_params = zip(periods, alphas)

#   more parameters for the fitting (like how many points we want to test for each, and what range of alpha and period we want it to go through)
num_alpha = 11
num_periods = 11
a_range = 0.000001 #+/- range of alpha to test
p_range = 0.25 #+/- range of periods to test

num_phase = 201 # this is just the number of phase points for the plot, won't really affect the run time
phase_points = np.linspace(0, 1.0, num_phase)
template_mags = ft.cosine_series_phase(coeff, phase_points, order)

fit_periods = []
fit_alphas = []
fit_chisq = []

for i in start_params:
    #   for each of the starting parameters, generate a cepheid light curve
    period = i[0]
    alpha = i[1]
    shifted_mags = ft.evolving_cosine(coeff, good_mjd, order, period, alpha, ft.linear_per_change, min_date=min(good_mjd))
    shifted_phases = good_mjd/period - np.fix(good_mjd/period)
    test_phases = np.asarray(pf.ceph_fold(good_mjd, period, alpha))

    #    add in some Gaussian noise
    noise = np.random.normal(0, 0.03, len(shifted_mags))
    mags_with_noise = shifted_mags + noise
    test_errors = np.ones(len(good_mag_v))*np.std(noise)

    #   make a plot to check that it's behaving correctly (can take this step out later)
    #plt.clf()
    #plt.scatter(shifted_phases, mags_with_noise, color='grey', alpha=0.3)
    #plt.scatter(test_phases, shifted_mags, color='red', alpha=0.3)
    #plt.gca().invert_yaxis()
    #plt.xlabel("Phase")
    #plt.ylabel("V-Band Magnitude")
    #plt.savefig("output/{}_fake_data.png".format(cepheid))
    #plt.show()

    best_alpha, best_period, chi_squares = pf.sec_change_per(good_mjd, mags_with_noise, coeff, period, test_errors, order, num_periods=num_periods, num_alpha=num_alpha)
    best_chi_square = pf.best_chi(chi_squares)
    best_phases = np.asarray(pf.ceph_fold(good_mjd, best_period, best_alpha))
    temp_mags = ft.cosine_series_phase(coeff, best_phases, order)
    temp_phases = np.asarray(pf.ceph_fold(good_mjd, best_period, best_alpha))
    
    #   for making the contour plot. probably not going to do this every time? or should i?
    alpha_range = np.linspace(-a_range, a_range, num_alpha)
    period_range = np.linspace(-p_range + period, p_range + period, num_periods)
    xv, yv = np.meshgrid(alpha_range, period_range)

    #plt.clf()
    #CS = plt.contour(xv, yv, chi_squares)
    #plt.clabel(CS, inline=1, fontsize=10)
    #plt.ylabel("Period")
    #plt.xlabel("Alpha")
    #plt.text(-0.00000059, 200, "Best Fit Chi-Square = {}".format(np.amin(chi_squares)), fontsize=14)
    #plt.savefig("output/{}_fake_chi_squares_contour.png".format(cepheid))
    #plt.show()

    #   repeat once for the plotting (or not)
    x, y = best_phases, mags_with_noise
    xor, yor = temp_phases, temp_mags

    #plt.clf()
    #plt.scatter(x, y, color='grey', alpha=0.3)
    #plt.scatter(xor, yor, color='red', alpha=0.3)
    #plt.gca().invert_yaxis()
    #plt.ylabel("V-Band Magnitudes")
    #plt.xlabel("Phase")
    #plt.savefig("output/{}_fake_final_pa.png".format(cepheid))
    #plt.show()

    fit_periods.append(best_period)
    fit_alphas.append(best_alpha)
    fit_chisq.append(best_chi_square)

#   output the results of the monte carlo to a .dat file
filename = 'output/mc_test_{}.dat'.format(n)

