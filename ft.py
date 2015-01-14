#   chuang
#   1.07.15
#
#   This module is supposed to contain everything useful. However, this is likely just the first iteration
#
####

#
#   First thing to do...write a function that gets the data for a certain Cepheid from the name
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.io.ascii as ascii
import pyfits
from PyAstronomy.pyTiming import pyPDM
from scipy import optimize
import heapq
import os.path

#
#   Gets ALL the data for that Cepheid (combining old ASAS, new ASAS, and IOMC)
#

def get_ceph_data(name):
    Cepheid_visual = np.genfromtxt('../iomc/Cepheid_visual.txt', dtype=None, names=['name', 'id', 'cr1', 'cr2', 'ra1', 'ra2', 'ra3', 'dec1', 'dec2', 'dec3', 'radec', 'decdec', 'cr11', 'type', 'cr22', 'cr33', 'm1', 'err', 'mmax', 'mmin', 'per', 'gal_l', 'gal_b', 'Schafly Av', 'SFD Av', 'mag/kpc'], skip_header=2, filling_values=-1.0)
    ids = Cepheid_visual['id']
    names = Cepheid_visual['name']
    id_numbers = np.arange(len(ids))
    cepheids = {}
    for i in id_numbers:
        cepheids[names[i]] = i
    per = Cepheid_visual['per']

    number = cepheids[name]
    period = per[number]
    star_data = Cepheid_visual[number]
    name = star_data[0]

    IOMC_dat = get_IOMC_data(number, period, ids[number])
    ASAS_dat = get_ASAS_data(name, period)

    if ASAS_dat.size != 1:
        print "ASAS data included"
        Ceph_data = np.concatenate([IOMC_dat, ASAS_dat], 0)
    else:
        print "no ASAS data"
        Ceph_data = IOMC_dat
    return Ceph_data

#
#   Gets all the ASAS data.
#
def get_ASAS_data(name, period):
    old_asas_file = "../asas_old/{}.txt".format(name)
    cepheid_filename = '_'.join([name[:-3],name[-3:]])
    new_asas_file = "../new_asas/{}.txt".format(cepheid_filename)
    if (os.path.isfile(new_asas_file) != False) & (os.path.isfile(old_asas_file) != False):
        oldasas = np.genfromtxt(old_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
        oldasas_date = oldasas['asasdate'] + 50000.0
        oldasas_epochs = np.fix(oldasas_date/period)
        oldasas_phase = oldasas_date/period - np.fix(oldasas_date/period)
        oldasas_grade = oldasas['grade']
        old_arr2 = np.arange(len(oldasas))
        old_asasgd = old_arr2[oldasas_grade[old_arr2] == 'A']
        old_len = len(old_asasgd)
        m2_old = oldasas['m2']
        m2_asasiomc_old = asas_old['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
        old_shape = (old_len, 1)
        old_phases = oldasas_phase[old_asasgd]
        old_m2 = m2_old[old_asasgd]
        dates = oldasas_date[old_asasgd]
        epochs = oldasas_epochs[old_asasgd]
    
        asas_new = np.genfromtxt(new_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
        new_asasdate = asas_new['asasdate'] + 50000.0
        new_asasepoch = np.fix(new_asasdate/period)
        new_asasphase = new_asasdate/period - np.fix(new_asasdate/period)
        grade_new = asas_new['grade']
        new_arr2 =  np.arange(len(asas_new))
        new_asasgd = new_arr2[grade_new[new_arr2] == 'A']
        new_len = len(new_asasgd)
        m2_new = asas_new['m2']
        m2_asasiomc_new = asas_new['m2'] - 0.03
        new_m2 = m2_new[new_asasgd] + 0.01
        n_m2 = m2_new[new_asasgd]
        new_phases = new_asasphase[new_asasgd]
        new_dates = new_asasdate[new_asasgd]
        new_epochs = new_asasepoch[new_asasgd]
        
        total_len = old_len + new_len
        arrshape = (total_len, 1)
    
        dates = np.concatenate([dates,new_dates],1)
        dates.shape = arrshape
        epochs = np.concatenate([epochs, new_epochs],1)
        epochs.shape = arrshape
        m2 = np.concatenate([old_m2, n_m2],1)
        m2.shape = arrshape
        phases = np.concatenate([old_phases, new_phases],1)
        phases.shape = arrshape
        
        print "new ASAS data available"
    
        ceph_arr = np.hstack([dates, m2, epochs, phases])
        return ceph_arr
    elif os.path.isfile(old_asas_file) != False:
        oldasas = np.genfromtxt(old_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
        oldasas_date = oldasas['asasdate'] + 50000.0
        oldasas_epochs = np.fix(oldasas_date/period)
        oldasas_phase = oldasas_date/period - np.fix(oldasas_date/period)
        oldasas_grade = oldasas['grade']
        old_arr2 = np.arange(len(oldasas))
        old_asasgd = old_arr2[oldasas_grade[old_arr2] == 'A']
        old_len = len(old_asasgd)
        arrshape = (old_len,1)
        m2_old = oldasas['m2']
        m2_asasiomc_old = oldasas['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
        old_phases = oldasas_phase[old_asasgd]
        old_phases.shape = arrshape
        old_m2 = m2_old[old_asasgd]
        old_m2.shape = arrshape
        dates = oldasas_date[old_asasgd]
        dates.shape = arrshape
        epochs = oldasas_epochs[old_asasgd]
        epochs.shape = arrshape
        
        print "no new ASAS data available"

        ceph_arr = np.hstack([dates, old_m2, epochs, old_phases])
        return ceph_arr
    else:
        return np.empty(1)
#
#   Gets the IOMC data for the Cepheid.
#   Final columns go in order of: mjd, mag_v, epoch, phase
#

def get_IOMC_data(number, period, id):
    btime = 51544.5
    filename = "../iomc/IOMC_{}.fits".format(id)
    hdulist = pyfits.open(filename)
    ceph_data = hdulist[1].data
    columns = hdulist[1].columns

    mag_v = ceph_data.field('mag_v')
    mjd = ceph_data.field('barytime') + btime
    errmag_v = ceph_data.field('errmag_v')
    exptime = ceph_data.field('exposure')
    pp = ceph_data.field('problems')

    newjd = mjd - min(mjd)

    phase = mjd/period - np.fix(mjd/period)
    sortidxs = np.argsort(phase)
    sorted_phase = phase[sortidxs]
    sorted_mag_v = mag_v[sortidxs]
    sorted_mjd = mjd[sortidxs]
    arr = np.arange(len(ceph_data))
    gd = arr[(exptime[arr] < 200) & (pp[arr] == np.median(pp))]
    length = len(gd)
    arrshape = (length, 1)
    mag_v = mag_v[gd]
    mag_v.shape = arrshape
    mjd = mjd[gd]
    mjd.shape = arrshape
    phase = phase[gd]
    phase.shape = arrshape
    epochs = np.fix(mjd/period)
    epochs.shape = arrshape
    
    ceph_arr = np.hstack([mjd, mag_v, epochs, phase])

    return ceph_arr

#
#   Helps find the initial parameters for the Fourier fit
#
def approximate_median(t, y):
    median = np.median(y)
    return median

def approximate_amplitude(t, y):
    amplitude = (max(y) - min(y))/2.
    return amplitude

#
#   Fold the data and then spit out the number of cycles needed, generally set to 2. Can give this the mjd and mag from the other parts and it will fold it. Then extract median will extract the median of the points and finally that can go into the fourier series.
#

def fold_data(t, y, period, num=2):


#
#   Extracts the median points from a set of data with lots of overlap
#

def extract_median(t, y):
    binned_t = []
    binned_mag = []
    t_diff = np.diff(t)
    end_ind = np.concatenate([np.where(t_diff > 0.04)[0] + 1, [len(t)]], 1)
    begin_ind = np.concatenate([[0], np.where(t_diff > 0.04)[0] + 1],1)
    bin_number = np.arange(0, len(end_ind))
    for i in bin_number:
        bin = t[begin_ind[i]:end_ind[i]]
        bin_mag = y[begin_ind[i]:end_ind[i]]
        t_median = np.median(bin)
        m_median = np.median(bin_mag)
        binned_t.append(t_median)
        binned_mag.append(m_median)
    
    binned_t = np.asarray(binned_t)
    binned_mag = np.asarray(binned_mag)
    return binned_t, binned_mag

#
#   Gets the Cepheid's period from the txt file
#

def get_period(name):
    Cepheid_visual = np.genfromtxt('../iomc/Cepheid_visual.txt', dtype=None, names=['name', 'id', 'cr1', 'cr2', 'ra1', 'ra2', 'ra3', 'dec1', 'dec2', 'dec3', 'radec', 'decdec', 'cr11', 'type', 'cr22', 'cr33', 'm1', 'err', 'mmax', 'mmin', 'per', 'gal_l', 'gal_b', 'Schafly Av', 'SFD Av', 'mag/kpc'], skip_header=2, filling_values=-1.0)
    ids = Cepheid_visual['id']
    names = Cepheid_visual['name']
    id_numbers = np.arange(len(ids))
    cepheids = {}
    for i in id_numbers:
        cepheids[names[i]] = i
    per = Cepheid_visual['per']

    number = cepheids[name]
    period = per[number]
    return period

#
#   Fits data with a Fourier Series. Before this works properly, there needs to be separate code to make the data easy to use for fitting...
#   Currently, it has issues fitting the phase (it has a ridiculously large, unjustified phase shift)
#   Basic idea is to try with one Fourier mode, and then use an f-test, add in another mode, another f-test, and stop adding in more modes when it turns out that the fit is not getting better
#   If you give the fourier fit an order, it will do a fourier fit to that order. Otherwise order is set to zero, meaning that it will keep adding in more fits until its not doing better.
#

def fourier_fit(mjd, mag, name, order=0):
    period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag) #this guess usually turns out to be really off
    zeroed_mag = mag - median # so that the data is more centered around zero...
    if order == 0:
        # fitfunc = lambda p, x: p[0]*cos(2*pi/p[1]*x+p[2]) + p[3]*x
        fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x + p[2]) + p[3]
        errfunc = lambda p, x, mag: fitfunc(p, x) - mag
        p0 = [amplitude, period, 0., 0.] # Initial guess, you'll probably have to write this into it somehow. Maybe have it call on pdm analysis first to get an approximate frequency? And something similar to get the median, etc. Start off by assuming a phase-shift of 0.
        p1, success = optimize.leastsq(errfunc, p0, args=(mjd,mag))
    
        time = np.linspace(mjd.min(), mjd.max(), 300)
        #time = np.linspace(53500, 54000, 300)
        plt.plot(mjd, mag, "ro", time, fitfunc(p1, time), "r-")
        plt.ylim([min(y), heapq.nlargest(2, y)[1]])
        #plt.xlim([53500, 54000])
        print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[0], '\n' + 'Period: ', p0[1], '\n'  + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[3], '\n'
        print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[0], '\n' + 'Period: ', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[3], '\n'
        
        plt.show()
    
        #   try higher orders
        
    else:
         print "still working on it"
    #   return fitfunc

#
#   Fit with a cosine
#

def cosine_func(params, xdata):
    return params[0]*np.cos(2*np.pi/params[1]*xdata + params[2]) + params[3]

#
#   Higher orders...try fixing the frequency first based off of the original cosine fit
#   To avoid confusion, instead of reordering the parameters, I will instead make the second element of parameter = 1, so that
#

def cosine_orders(params, xdata, order, period):
    return params[0]*np.cos(2*np.pi/order*period*xdata*params[1] + params[2]) + params[3]

#
#   Error function
#

def error_func(params, function, xdata, ydata):
    return function(params, xdata) - ydata

#
#   Another attempt at an error function
#


def error_series(params, xdata, ydata, order, period):
    return cosine_series(params, xdata, order, period) - ydata

#
#   New approach. The input parameters may just be the guess for the original mode
#

def cosine_series(params, xdata, order, period):
    fund_mode = params[0]*np.cos(2*np.pi/period*xdata + params[1]) + params[2]
    if order > 1:
        sum_func = fund_mode
        for i in (np.arange(order-1) + 2):
            sum_func = fund_mode + params[i+1]*np.cos(2*np.pi/(2*period)*xdata + params[i+2])
    else:
        sum_func = fund_mode
    return sum_func

#
#   New version of fourier fit. Currently is able to fit to a certain mode that you input. Need to change it to so that it will fit to any mode, depending on the errors.
#

def fourier_series(mjd, mag, name, order=1.):
    #   initial guess for period
    period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    p0 = [amplitude, 0., median]
    
    p1, success = optimize.leastsq(error_series, p0, args=(mjd, mag, order, period))
    errors = np.sum(np.power(error_series(p1, mjd, mag, order, period), 2))
    add_order = True
    
    while (add_order == True) & (order < 6):
        order = order + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(error_series, p00, args=(mjd, mag, order, period))
        errors2 = np.sum(np.power(error_series(p2, mjd, mag, order, period), 2))
        if errors2 < errors:
            p1 = p00*1.
            errors = errors2*1.
        else:
            add_order = False
            order = order - 1
    
    # time = np.linspace(mjd.min(), mjd.max(), 300)
    time = np.linspace(53500, 53600, 300)
    plt.plot(mjd, mag, "ro", time, cosine_series(p1, time, order, period), "r-")
    plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[0],'\n' + 'Phase Shift: ', p0[1], '\n' + 'Offset: ', p0[2], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[0], '\n'  + 'Phase Shift: ', p1[1], '\n' + 'Offset: ', p1[2], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', order
    print 'PERIOD: ', period
    print p1
    print p2

    plt.show()

#fixed_per = p1[1]
#p11 = [p1[0], 1., p1[2], p1[3]]
#print p1
#order = 2
#function2 = cosine_func(p1, mjd)
#p2, success = optimize.leastsq(error_func2, p11, args=(cosine_orders, function2, mjd, mag, order, fixed_per))
#errors2 = np.sum(np.power(error_func2(cosine_orders, function2, p2, mjd, mag, order, period), 2))
