#	chuang
#	12.10.14
#
#	module to plot cepheids in a variety of ways -- give it an array of times, magnitudes, and a period and it will plot in "standard" fashion with two periods shown

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#   Makes a plot in "the usual" way where the periods are doubled
#
#   Parameters:
#   t = the dates
#   y = the magnitudes
#   period = the cepheid's period
#   name = the cepeheid's name
#   t,y are numpy arrays, and period is a float

def standard_plot(t, y, period, name='none'):

    phase = (t/period)*1.0 - np.fix(t/period)
    sortidxs = np.argsort(phase)
    sorted_phase = phase[sortidxs]
    sorted_y = y[sortidxs]
    sorted_dates = t[sortidxs]
    pp = np.concatenate([sorted_phase, sorted_phase+1.0], 1)
    yy = np.concatenate([sorted_y, sorted_y], 1)
    mi = min(y)
    plt.clf()
    plt.scatter(pp, yy, alpha=0.7, color='grey')
    plt.gca().invert_yaxis()
    plt.text(-0.4, mi, "period = {} days".format(period), fontsize=14)
    plt.xlabel("phase")
    plt.ylabel("V-band magnitude")
    if name !='none':
        plt.title("{}".format(name))
        plt.savefig("../output/{}_standard.png".format(name))
    plt.show()

#   Returns the array to let you make a standard plot, but doesn't actually make the plot (which standard_plot does)
#
#   Parameters:
#   t = the dates
#   y = the magnitudes
#   period = the cepheid's period

def splot_array(t, y, period):
    
    phase = (t/period)*1.0 - np.fix(t/period)
    sortidxs = np.argsort(phase)
    sorted_phase = phase[sortidxs]
    sorted_y = y[sortidxs]
    sorted_dates = t[sortidxs]
    pp = np.concatenate([sorted_phase, sorted_phase+1.0], 1)
    yy = np.concatenate([sorted_y, sorted_y], 1)
    return pp, yy

#   Makes a plot in which the epochs are color-coded using heatmap coloring (earliest in red, latest in blue)
#
#   Parameters:
#   t = the dates
#   y = magnitudes
#   period = the cepheid's period

def epoch_plot(t, y, period, name='none'):

    phase = (t/period)*1.0 - np.fix(t/period)
    sortidxs = np.argsort(phase)
    sorted_phase = phase[sortidxs]
    sorted_y = y[sortidxs]
    sorted_dates = t[sortidxs]
    pp = np.concatenate([sorted_phase, sorted_phase+1.0], 1)
    yy = np.concatenate([sorted_y, sorted_y], 1)
    epochs = np.fix(sorted_dates/period)
    epoch_list = set(epochs)
    mi = min(y)
    plt.clf()
    x = 0
    for i in epoch_list:
        period_date_idxs = np.asarray(np.where(epochs == i))
        period_mags = sorted_y[period_date_idxs]
        period_dates = sorted_dates[period_date_idxs]
        period_phases = sorted_phase[period_date_idxs]
        pp = np.concatenate([period_phases, period_phases + 1.0], 1)
        yy = np.concatenate([period_mags, period_mags], 1)
        plt.scatter(pp, yy, color=cm.jet(1.*x/len(epoch_list)), alpha=0.7, edgecolors='none')
        x +=1
    plt.gca().invert_yaxis()
    plt.text(-0.4, mi, "period = {} days".format(period), fontsize=14)
    plt.xlabel("phase")
    plt.ylabel("V-band magnitude")
    if name !='none':
        plt.title("{}".format(name))
        plt.savefig("../output/{}_epoch.png".format(name))
    plt.show()


