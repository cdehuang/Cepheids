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

folded_phase, mags = ft.fold_data(mjd, mag_v, period)
ph_med, run_med = ft.get_median(folded_phase, mags)

coeff, ord = ft.fourier_series(ph_med*period, run_med, name)
#   ft.fourier_series(folded_phase, mags, name) <-- doesn't really work

#   Next task is to estimate the error. This will be done by the method that Richard suggested yesterday -- moving the
