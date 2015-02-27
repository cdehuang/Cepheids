#
#	chuang
#	02.24.2015
#
#

"""
	redoing the fit with just one period input and using weighted least squares fit.
"""

import numpy as np
import matplotlib.pyplot as plt

import ft
import cepheid_plotting

cepheid = 'xycar'

#   get the data
ceph_data = ft.get_sorted_data(cepheid)
mjd = ceph_data[:,0]
mag_v = ceph_data[:,1]
epochs = ceph_data[:,2]
phases = ceph_data[:,3]
mag_v_err = ceph_data[:,4]

#coeff, order = ft.fourier_series(mjd, mag_v, cepheid)
coeff, order = ft.fs_wtlsq(mjd, mag_v, cepheid, mag_v_err, order=7)

#   plot with the fit
period = ft.get_period(cepheid)
phs, ydata, errs = ft.fold_data(mjd, mag_v, period, mag_v_err)

dates = np.linspace(0, 2*period, 1000)

ft_y = ft.cosine_series(coeff, dates, 7, period)

plt.scatter(phs, ydata, color='grey', alpha=0.3)
plt.plot(dates/period, ft_y, color='red', alpha=0.7, linewidth=3)
plt.gca().invert_yaxis()
plt.title("{}".format(cepheid))
plt.ylabel("V-band Magnitude")
plt.xlabel("Phase")
plt.savefig("../output/{}_with_fourier_fit".format(cepheid))
plt.show()



if __name__=="__main__":
    print "Nothing here yet"
