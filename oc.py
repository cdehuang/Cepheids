#	chuang
#	1.19.15
#
#	O-C Diagram program. Does a Fourier fit to one (folded) period of the data and then it tries to calculate what part of the phase that it's at using just the magnitude. Then try to find what part of the period you're looking at after doing this (by picking the period that is closer to your data points...note that this method gets pretty tricky if you happen to be at one of the turnaround points...)

import numpy as np
import ft
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
    # return array[idx]

def berdnikov(M0, period, q, epoch):
    berdy = M0 + period*epoch + q*(epoch**2)
    return berdy

def findintersect(fun1, fun2, x0):
    return fsolve(lambda x: fun1(x) - fun2(x), x0)

def hline(xval, yval):
    return yval

#   should probably try to see how this compares with our data

name = 'svul' # <-- this should be the sole input parameter for now. later on the other input parameter will be the date of observation
cepheid_data = ft.get_ceph_data(name)
mjd = cepheid_data[:,0]
mag_v = cepheid_data[:,1]
epoch = cepheid_data[:,2]
phase = cepheid_data[:,3]
err = cepheid_data[:,4]
period = ft.get_period(name)

#   Try sorting everything first, so now all the data has been sorted so that it's in chronological order

sort_inds = np.argsort(mjd)
mjd = mjd[sort_inds]
mag_v = mag_v[sort_inds]
epoch = epoch[sort_inds]
phase = phase[sort_inds]
err = err[sort_inds]

#   note that a seventh-order polynomial fit was used for this just so that I could get a really good idea of the shape...

folded_phase, mags, errors = ft.fold_data(mjd, mag_v, period, err=err)
ph_med, run_med, err_med = ft.get_median(folded_phase, mags, errors=errors)
coeff, ord = ft.fourier_series(ph_med*period, run_med, name, order=3)
phase_shifts = np.linspace(-0.5, 0.5, 100)

#   fit function
order = len(coeff)/2
xval = np.linspace(0, 1., 201)
yval = ft.cosine_series(coeff, xval*period, order, period)

#   need to toss out areas where it is contentious maybe? Or at least toss out ones where the magnitude is much much greater (or less than) it's supposed to be. You could also do better if you zeroed the phases on the lowest point (in terms of matching them up). Calculating the slope might be difficult with this number of points

phase_points = []
phase_idxs = []
discarded = []

#   WARNING: Has been modified so that it only works if everything is sorted
#   Maybe a better idea is to draw a line through and find where it intersects? And then try to pick the phase that matches the best?
"""

for i in range(len(mag_v)):
    m = mag_v[i]
    p = phase[i]
    result = findintersect(hline, )

"""

for i in range(len(mag_v)):
    m = mag_v[i]
    p = phase[i]
    id = find_nearest(m, yval)
    if np.abs(yval[id]-m) > 0.05:
        print "probably bad data. magnitude: ", m
        discarded.append(i)
    elif np.abs(xval[id] - p) > 0.2:
        if ((mjd[i+1] - mjd[i]) < period) & (i+1 not in discarded):
            phase_diff = (mjd[i+1] - mjd[i])/period
            print "hold that thought"
        
        elif ((mjd[i] - mjd[i-1]) < period) & (i-1 not in discarded):
            print "close data available: ",
        
        else:
            print "bad fit. data index: ", id
            print "fit index: ", i
            discarded += 1
    else:
        phase_points.append(xval[id])
        phase_idxs.append(i)

print "discarded", discarded, "data points"
print len(phase_points), "good data points"


#   here are the mjd for the phase measurements that matched up
good_mjd = mjd[phase_idxs]
good_epochs = epoch[phase_idxs]

#   now we have to calculate what mjd I would actually have if I multiplied the phase by the epoch

calc_mjd = [period*(good_epochs[i] + phase_points[i]) for i in range(len(good_epochs))]

diff = calc_mjd - good_mjd
plt.scatter(good_mjd, diff)
plt.savefig("{}_oc.png".format(name))
plt.show()

