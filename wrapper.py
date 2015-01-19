#   chuang
#   1.15.15
#
#   I expect that this will probably undergo some more changes, hence the kind of dumb name.
#
#
#   Basically, this code's job is to just run the fourier transform stuff (which is still under improvement at the moment) and calculate the uncertainties. Though I'll probably have some separate module for the uncertainty part? Haven't really decided on that yet though.
#
#

import numpy as np
import matplotlib.pyplot as plt
import ft
from scipy import stats

#   Calculate the residuals

def residuals(observed, expected):
    return observed - expected

#   Calculate chi-squared --> the errors need to be assessed in a different way, but for now just get the errors on those specific data points I suppose?

def chi_square(resid, errors):
    cs = np.sum(np.power(resid,2.)/np.power(errors,2.))
    return cs

name = 'svul' # <-- this should be the sole input parameter for now. later on the other input parameter will be the date of observation
cepheid_data = ft.get_ceph_data(name)
mjd = cepheid_data[:,0]
mag_v = cepheid_data[:,1]
epoch = cepheid_data[:,2]
phase = cepheid_data[:,3]
err = cepheid_data[:,4]
period = ft.get_period(name)

#
#   Not sure yet if the fitting works better with one period or with multiple periods. So in that case, I don't know if I would want to use the folded data into the fourier fitting thing or just the regular data. For now we will assume that having multiple periods means that it fits better...probably doesn't really though
#

folded_phase, mags, errors = ft.fold_data(mjd, mag_v, period, err=err)
ph_med, run_med, err_med = ft.get_median(folded_phase, mags, errors=errors)

coeff, ord = ft.fourier_series(ph_med*period, run_med, name)
#   ft.fourier_series(folded_phase, mags, name) <-- doesn't really work

#   Next task is to estimate the error. This will be done by the method that Richard suggested yesterday -- moving the light curve around from the average by a phase and then calculating what the chi-squared is at each location.
#
#   Error calculation

phase_shifts = np.linspace(-0.5, 0.5, 100)
chi_sq_array = []
for i in phase_shifts:
    shifted_phases = (ph_med+i)*period

    expected = ft.cosine_series(coeff, ph_med*period, ord, period)
    observed = ft.cosine_series(coeff, shifted_phases, ord, period)
    
    resid = residuals(observed, expected)

    output = chi_square(resid, err_med)
    chi_sq_array.append(output)

plt.scatter(phase_shifts, chi_sq_array)
plt.ylabel("$\chi^2$")
plt.xlabel("Phase Shift")
plt.title("{}".format(name))
#plt.savefig()
plt.show()



"""
#   Plot the shift
phase_shift = 0.1
plt.plot(ph_med*period, ft.cosine_series(coeff, ph_med*period, ord, period))
plt.plot(ph_med*period, ft.cosine_series(coeff, shifted_phases, ord, period), color='red')
plt.scatter(ph_med*period, run_med, alpha=0.4, color='grey')
plt.show()
"""
