#
#	chuang
#	03.06.2015
#

import numpy as np
import matplotlib.pyplot as plt

import phase_folding as pf
import ft
#phase_folding should import everything else, so it's probably ok after this

"""
    This is some example code that shows how the phase_folding module works.
    
    Basically, it makes some fake data with some linear period change and sees if my code can reproduce it.
    
    The program should be able to match everything 'exactly' in this case because the template is right, all I have done is add some noise.
    
    Should also do a monte Carlo where you start of with a bunch of different periods and alphas (maybe random uniform distribution)? And then you spit out the chi_sq, the period it finds, and the the alpha for each thing.
    
    This is, of course, assuming that I have a model that perfectly describes the data. If the model doesn't perfectly describe the data, then I have to run a separate set of trials for that (I think).
    
    Monte Carlo
"""

cepheid = 'xycar'

ceph_data = ft.get_sorted_data(cepheid)
mjd = ceph_data[:,0]
mag_v = ceph_data[:,1]
epochs = ceph_data[:,2]
phases = ceph_data[:,3]
mag_v_err = ceph_data[:,4]

#   do a crappy clipping of the really bad points--fix this later. need to write a sigma-clipping program in the pf module.
good_ind = pf.simple_sigma_clip(mag_v, phases)
good_ceph_dat = ceph_data[good_ind]
good_mjd = mjd[good_ind]
good_mag_v = mag_v[good_ind]
good_epochs = epochs[good_ind]
good_phases = phases[good_ind]
good_mag_v_err = mag_v_err[good_ind]

#   starting parameters of the fake data
period = ft.get_period(cepheid) + 0.001
alpha = 0.0000001

#   make a template cepheid--based off of a real cepheid
coeff, order = ft.fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=7)

#   create cepheid data with some complex phase that is based on the template cepheid.
shifted_mags = ft.evolving_cosine(coeff, good_mjd, order, period, alpha, ft.linear_per_change, min_date = min(mjd))
shifted_phases = good_mjd/period - np.fix(good_mjd/period)
test_phases = np.asarray(pf.ceph_fold(good_mjd, period, alpha))

#   add in some Gaussian noise
noise = np.random.normal(0, 0.03, len(shifted_mags))
mags_with_noise = shifted_mags + noise
test_errors = np.ones(len(good_mag_v))*np.std(noise)

#   plot to make sure that it's actually doing the right things
plt.clf()
plt.scatter(shifted_phases, mags_with_noise, color='grey', alpha=0.3)
plt.scatter(test_phases, shifted_mags, color='red', alpha=0.3)
plt.gca().invert_yaxis()
plt.xlabel("Phase")
plt.ylabel("V-Band Magnitude")
plt.savefig("output/{}_fake_data.png".format(cepheid))
plt.show()

#   now set up the fitting with the number of points to test
num_alpha = 21
num_periods = 21
a_range = 0.000001
p_range = 0.25
num_phase = 201 # just the phase_points for the plot
phase_points = np.linspace(0, 1.0, num_phase)
template_mags = ft.cosine_series_phase(coeff, phase_points, order)

#   find the change that works the best
best_alpha, best_fit_period, chi_squares = pf.sec_change_per(good_mjd, mags_with_noise, coeff, period, test_errors, order, num_periods=num_periods, num_alpha=num_alpha)
best_phases = np.asarray(pf.ceph_fold(good_mjd, best_fit_period, best_alpha))
temp_mags = ft.cosine_series_phase(coeff, best_phases, order)
temp_phases = np.asarray(pf.ceph_fold(good_mjd, best_fit_period, best_alpha))
alpha_range = np.linspace(-a_range, a_range, num_alpha)
period_range = np.linspace(-p_range + period, p_range + period, num_periods)
xv, yv = np.meshgrid(alpha_range, period_range)

plt.clf()
CS = plt.contour(xv, yv, chi_squares)
plt.clabel(CS, inline=1, fontsize=10)
#plt.xlim([-0.000006, 0.000006])
plt.ylabel("Period")
plt.xlabel("Alpha")
plt.text(-0.00000059, 200, "Best Fit Chi-Square = {}".format(np.amin(chi_squares)), fontsize=14)
plt.savefig("output/{}_fake_chi_squares_contour.png".format(cepheid))
plt.show()

#   repeat once for the plotting (or not)
x, y = best_phases, mags_with_noise
xor, yor = temp_phases, temp_mags
#xor, yor = pf.duplicate(shifted_phases, mags_with_noise)

plt.clf()
plt.scatter(x, y, color='grey', alpha=0.3)
plt.scatter(xor, yor, color='red', alpha=0.3)
#plt.scatter(best_phases, good_mag_v, color='grey', alpha=0.3)
#plt.scatter(good_phases, good_mag_v, color='red', alpha=0.3)
plt.gca().invert_yaxis()
plt.ylabel("V-Band Magnitudes")
plt.xlabel("Phase")
plt.savefig("output/{}_fake_final_pa.png".format(cepheid))
plt.show()
