#   chuang
#   1.07.15
#

"""
 This module is supposed to contain everything useful for the Cepheid data analysis. However, this is likely just the first iteration
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.io.ascii as ascii
import pyfits
from PyAstronomy.pyTiming import pyPDM
from scipy import optimize
import heapq
import os.path
import math

def get_ceph_data(name, mag_shift=0.0):
    """
        Accidentally deleted the original get ceph data, so here it is again.
    """
    
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
    ASAS_dat = get_ASAS_data(name, period, mag_shift=mag_shift)
    
    #print ASAS_dat.shape
    #print IOMC_dat.shape
    
    if ASAS_dat.size != 1:
        print "ASAS data included"
        Ceph_data = np.concatenate([IOMC_dat, ASAS_dat], 0)
    else:
        print "no ASAS data"
        Ceph_data = IOMC_dat

    return Ceph_data

def make_data_file(name, fname, mag_shift=0.0):
    """
        Does the same thing as get_ceph_data, except that instead of returning a light curve it prints the output to a dat file instead.
        
        Returns the filename of the output file.
        """
    
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
    ASAS_dat = get_ASAS_data(name, period, mag_shift=mag_shift)
    
    #print ASAS_dat.shape
    #print IOMC_dat.shape
    
    if ASAS_dat.size != 1:
        print "ASAS data included"
        Ceph_data = np.concatenate([IOMC_dat, ASAS_dat], 0)
    else:
        print "no ASAS data"
        Ceph_data = IOMC_dat
    
    #   Sort all of the data into chronological order
    sortidxs = np.argsort(Ceph_data[:,0])
    Ceph_data = Ceph_data[sortidxs]
    
    #   convert into a record array
    Ceph_rec = np.core.records.fromarrays(Ceph_data.transpose(), names='mjd, V-band magnitude, epochs, phases, errors')
    file_name = 'output/' + fname + '_{}.dat'.format(name)
    ascii.write(Ceph_rec, file_name)
    
    return file_name

def make_ceph_lc():
    """
        Makes a Cepheid light curve (so it only keeps the magnitude, error, and days). But it also allows to pick which band you want. The other code should be adjusted for this too at some point.
    """
    return "TODO"

def get_sorted_data(name, mag_shift=0.0):
    """
    Gets all of the available data for that Cepheid (new ASAS, old ASAS, and IOMC)
    """
    
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
    ASAS_dat = get_ASAS_data(name, period, mag_shift=mag_shift)

    print ASAS_dat.shape
    print IOMC_dat.shape

    if ASAS_dat.size != 1:
        print "ASAS data included"
        Ceph_data = np.concatenate([IOMC_dat, ASAS_dat], 0)
    else:
        print "no ASAS data"
        Ceph_data = IOMC_dat

    sortidxs = np.argsort(Ceph_data[:,0])
    Ceph_data = Ceph_data[sortidxs]
    return Ceph_data


def get_ASAS_data(name, period, mag_shift=0.0):
    """
    Gets all of the ASAS data available for the Cepheid
    """
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
        oldasas_e2 = oldasas['e2']
        m2_asasiomc_old = oldasas['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
        old_shape = (old_len, 1)
        old_phases = oldasas_phase[old_asasgd]
        old_m2 = m2_old[old_asasgd] - mag_shift
        dates = oldasas_date[old_asasgd]
        epochs = oldasas_epochs[old_asasgd]
        e2 = oldasas_e2[old_asasgd]
    
        asas_new = np.genfromtxt(new_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
        new_asasdate = asas_new['asasdate'] + 50000.0
        new_asasepoch = np.fix(new_asasdate/period)
        new_asasphase = new_asasdate/period - np.fix(new_asasdate/period)
        grade_new = asas_new['grade']
        new_arr2 =  np.arange(len(asas_new))
        new_asasgd = new_arr2[grade_new[new_arr2] == 'A']
        new_len = len(new_asasgd)
        m2_new = asas_new['m2']
        e2_new = asas_new['e2']
        m2_asasiomc_new = asas_new['m2'] - 0.03
        new_m2 = m2_new[new_asasgd] + 0.01
        n_m2 = m2_new[new_asasgd] - mag_shift
        new_phases = new_asasphase[new_asasgd]
        new_dates = new_asasdate[new_asasgd]
        new_epochs = new_asasepoch[new_asasgd]
        new_e2 = e2_new[new_asasgd]
        
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
        e2 = np.concatenate([e2, new_e2], 1)
        e2.shape = arrshape
        
        print "new ASAS data available"
    
        ceph_arr = np.hstack([dates, m2, epochs, phases, e2])
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
        oldasas_e2 = oldasas['e2']
        m2_asasiomc_old = oldasas['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
        old_phases = oldasas_phase[old_asasgd]
        old_phases.shape = arrshape
        old_m2 = m2_old[old_asasgd] - mag_shift
        old_m2.shape = arrshape
        dates = oldasas_date[old_asasgd]
        dates.shape = arrshape
        epochs = oldasas_epochs[old_asasgd]
        epochs.shape = arrshape
        e2 = oldasas_e2[old_asasgd]
        e2.shape = arrshape
        
        print "no new ASAS data available"

        ceph_arr = np.hstack([dates, old_m2, epochs, old_phases, e2])
        return ceph_arr
    else:
        return np.empty(1)

#
#   Gets only the new ASAS data
#

def get_new_ASAS_data(name, period, mag_shift=0.0):
    """
    Gets only the new ASAS data for the Cepheid
    
    INPUTS
    name: Name of the cepheid with no spaces (string)
    period: the Cepheid's period (float)
    mag_shift: the magnitude shift between the ASAS data and whatever other dataset (generally the IOMC data) (float)
    
    OUTPUTS
    An array containing the cepheid data. First column of the array are the Julian dates, the second column is the magnitudes, third the epochs, fourth the phases, and the last column the error
    """
    cepheid_filename = '_'.join([name[:-3],name[-3:]])
    new_asas_file = "../new_asas/{}.txt".format(cepheid_filename)
    if (os.path.isfile(new_asas_file) != False):
        asas_new = np.genfromtxt(new_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
        new_asasdate = asas_new['asasdate'] + 50000.0
        new_asasepoch = np.fix(new_asasdate/period)
        new_asasphase = new_asasdate/period - np.fix(new_asasdate/period)
        grade_new = asas_new['grade']
        new_arr2 =  np.arange(len(asas_new))
        new_asasgd = new_arr2[grade_new[new_arr2] == 'A']
        new_len = len(new_asasgd)
        m2_new = asas_new['m2']
        e2_new = asas_new['e2']
        m2_asasiomc_new = asas_new['m2'] - 0.03
        new_m2 = m2_new[new_asasgd] + 0.01
        n_m2 = m2_new[new_asasgd] - mag_shift
        new_phases = new_asasphase[new_asasgd]
        new_dates = new_asasdate[new_asasgd]
        new_epochs = new_asasepoch[new_asasgd]
        new_e2 = e2_new[new_asasgd]
        dates.shape = arrshape
        epochs.shape = arrshape
        m2.shape = arrshape
        phases.shape = arrshape
        e2.shape = arrshape

        ceph_arr = np.hstack([dates, new_m2, epochs, new_phases, e2])
        return ceph_arr


def get_old_ASAS_data(name, period, mag_shift=0.0):
    """
    Gets only the old ASAS data
    
    INPUTS
    name: Name of the cepheid with no spaces (string)
    period: the Cepheid's period (float)
    mag_shift: the magnitude shift between the ASAS data and whatever other dataset (generally the IOMC data) (float)
    
    OUTPUTS
    An array containing the cepheid data. First column of the array are the Julian dates, the second column is the magnitudes, third the epochs, fourth the phases, and the last column the error
    """
    old_asas_file = "../asas_old/{}.txt".format(name)
    if os.path.isfile(old_asas_file) != False:
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
        oldasas_e2 = oldasas['e2']
        m2_asasiomc_old = oldasas['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
        old_phases = oldasas_phase[old_asasgd]
        old_phases.shape = arrshape
        old_m2 = m2_old[old_asasgd] - mag_shift
        old_m2.shape = arrshape
        dates = oldasas_date[old_asasgd]
        dates.shape = arrshape
        epochs = oldasas_epochs[old_asasgd]
        epochs.shape = arrshape
        e2 = oldasas_e2[old_asasgd]
        e2.shape = arrshape
        ceph_arr = np.hstack([dates, old_m2, epochs, old_phases, e2])
        return ceph_arr

#
#   Gets the IOMC data for the Cepheid.
#   Final columns go in order of: mjd, mag_v, epoch, phase, err
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
    errmag_v = errmag_v[gd]
    errmag_v.shape = arrshape
    
    ceph_arr = np.hstack([mjd, mag_v, epochs, phase, errmag_v])

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

def fold_data(t, y, period, err=0, num=2):
    phase = t/period - np.fix(t/period)
    phases = phase*1.
    ydata = y*1.
    errors = err*1.
    print num
    x = 0
    if (num > 1) & (isinstance(err, int) != True):
        for i in (np.arange(num-1) + 2):
            x += 1
            print x
            add_phase = phase + (i-1.)
            phases = np.concatenate([phases, add_phase])
            ydata = np.concatenate([ydata, y])
            errors = np.concatenate([errors, err])
        return phases, ydata, errors
    elif (num > 1):
        for i in (np.arange(num-1) + 2):
            x += 1
            print x
            add_phase = phase + (i-1.)
            phases = np.concatenate([phases, add_phase])
            ydata = np.concatenate([ydata, y])
        return phases, ydata
    elif (isinstance(err, int) != True):
        return phases, ydata, errors
    else:
        return phases, ydata

# ipython notebook for the steve stuff
# Simpler median extraction--just bin by bins of set size. Unless the number of bins is specified, it will just do about 50 bins per cycle
#
#   Empty bins get taken out--maybe want to make a minimum bin population also? Not sure what the best way to do this is, if it's better to take the median over two bins or what
#   Get an approximate measurement of error by adding up the errors and dividing by the number of things in a bin. Fix this to make it actually sound later.

def get_median2(phases, ydata, errors=0, num_bins=50):
    #   total number of bins = bins
    num_phases = max(phases)
    nbins = int(num_bins*math.ceil(num_phases))
    
    #   bin data
    bins = np.linspace(0, math.ceil(num_phases), nbins)
    delta = bins[1] - bins[0]
    idx = np.digitize(phases, bins)
    bin_indices = range(nbins)
    for i in bin_indices:
        if len(ydata[idx == i]) > 5:
            running_median.append(np.median(ydata[idx ==i]))
            bin_center = bins[i] - delta/2
        elif (len(ydata[idx == i]) + len(ydata[idx == i + 1]) > 5):
            print "think more about this"
        else:
            print "think more about this too"
    #   running_median = np.asarray([np.median(ydata[idx == k]) for k in range(nbins)])
    phase_medians = bins - delta/2
    
    gd = np.where(np.logical_not(np.isnan(running_median)))[0]
    gd_pm = phase_medians[gd]
    gd_rm = running_median[gd]
    
    if (errors.size > 1):
        err = np.asarray([np.average(errors[idx == k]) for k in gd])
        
        #   necessary to slice before returning to remove the nan elements
        
        return gd_pm, gd_rm, err
    else:
        return gd_pm, gd_rm, 0
"""
    #   plot as a test
    plt.scatter(phases, ydata, color='k', alpha=0.2, s=2)
    plt.plot(bins-delta/2, running_median, 'r--', lw=4, alpha=.8)
    plt.axis('tight')
    plt.show()
"""

def get_median(phases, ydata, errors=0, num_bins=50):
    #   total number of bins = bins
    num_phases = max(phases)
    nbins = int(num_bins*math.ceil(num_phases))
    
    #   bin data
    bins = np.linspace(0, math.ceil(num_phases), nbins)
    delta = bins[1] - bins[0]
    idx = np.digitize(phases, bins)
    running_median = np.asarray([np.median(ydata[idx == k]) for k in range(nbins)])
    phase_medians = bins - delta/2
    
    gd = np.where(np.logical_not(np.isnan(running_median)))[0]
    gd_pm = phase_medians[gd]
    gd_rm = running_median[gd]
    
    if (errors.size > 1):
        err = np.asarray([np.average(errors[idx == k]) for k in gd])
        
        #   necessary to slice before returning to remove the nan elements
        
        return gd_pm, gd_rm, err
    else:
        return gd_pm, gd_rm, 0


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

def get_IOMC_id(name):
    """
        Gets the IOMC ID number of the Cepheid for when you want to use the get_IOMC_data() function
    """
    Cepheid_visual = np.genfromtxt('../iomc/Cepheid_visual.txt', dtype=None, names=['name', 'id', 'cr1', 'cr2', 'ra1', 'ra2', 'ra3', 'dec1', 'dec2', 'dec3', 'radec', 'decdec', 'cr11', 'type', 'cr22', 'cr33', 'm1', 'err', 'mmax', 'mmin', 'per', 'gal_l', 'gal_b', 'Schafly Av', 'SFD Av', 'mag/kpc'], skip_header=2, filling_values=-1.0)
    ids = Cepheid_visual['id']
    names = Cepheid_visual['name']
    id_numbers = np.arange(len(ids))
    cepheids = {}
    for i in id_numbers:
        cepheids[names[i]] = i
    per = Cepheid_visual['per']
    
    number = cepheids[name]
    id_number = ids[number]
    return id_number

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
        
        plt.ylabel("V-band Magnitude")
        plt.xlabel("Days")
        plt.title("{}".format(name))
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
#   Another attempt at an error function--with weights
#

def weighted_error_series(params, xdata, ydata, order, period, errors):
    return (cosine_series(params, xdata, order, period) - ydata)*(1/errors)

#
#  Switched around the order of the parameters so that they make more sense.
#

def cosine_series(params, xdata, order, period):
    sum_func = params[0]
    for i in np.arange(order) + 1:
        sum_func = sum_func + params[2*i-1]*np.cos(2*i*np.pi/(period)*xdata + params[2*i])
    return sum_func

def linear_per_change(period, alpha, mjd, min_date=0.0):
    return period + alpha*(mjd-min_date)

def evolving_cosine(params, xdata, order, start_period, alpha, per_func, min_date=0.0):
    #min_date should be min(mjd)
    sum_func = params[0]
    for i in np.arange(order) + 1:
        sum_func = sum_func + params[2*i-1]*np.cos(2*i*np.pi/(per_func(start_period, alpha, xdata, min_date))*xdata + params[2*i])
    return sum_func

def cosine_series_phase(params, phase, order):
    sum_func = params[0]
    for i in np.arange(order) + 1:
        sum_func = sum_func + params[2*i-1]*np.cos(2*i*np.pi*phase + params[2*i])
    return sum_func
#
#   Yet another attempt at an error function--with weights, but this time as 15.7.10 from Numerical Recipes
#


#def weighted_lorentz(params, xdata, ydata, order, period, errors):
    """
        Warning, not complete yet
    """
#    return (cosine_series(params, xdata, order, period) - ydata)*(1/errors))

#
#   Sigma clipping
#

def sigma_clip(data, sigma=3):
    
    return good_indices

"""
    if order > 1:
        sum_func =
        for i in (np.arange(order-1) + 2):
            sum_func = sum_func + params[i*2-1]*np.cos(2*np.pi/(2*period)*xdata + params[i*2])
    else:
        sum_func = params[1]*np.cos(2*np.pi/period*xdata + params[2]) + params[0]
"""


#
#   New version of fourier fit. Currently is able to fit to a certain mode that you input.
#   Update with BIC penalty. Chose 7 because that seems to be the fit order that most people use. Changed to 4 because it was clearly overfitting some of them.
#

def fourier_series(mjd, mag, name, order=7):
    #   initial guess for period
    ord = 1.
    period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    p0 = [median, amplitude, 0.]
    
    p1, success = optimize.leastsq(error_series, p0, args=(mjd, mag, ord, period))
    errors = np.sum(np.power(error_series(p1, mjd, mag, ord, period), 2))
    add_order = True
    
    while (add_order == True) & (ord < order):
        ord = ord + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(error_series, p00, args=(mjd, mag, ord, period))
        errors2 = np.sum(np.power(error_series(p2, mjd, mag, ord, period), 2))
        if errors2 < errors:
            p1 = p00*1.
            errors = errors2*1.
        else:
            add_order = False
#   order = order - 1
    
    time = np.linspace(mjd.min(), mjd.max(), 300)
    #   time = np.linspace(53500, 53600, 300)
    plt.scatter(mjd, mag, alpha=0.4, color='grey')
    plt.plot(time, cosine_series(p2, time, ord, period), "r-")
    plt.gca().invert_yaxis()
    #   plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    #   plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[1],'\n' + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[0], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[0], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', ord
    print 'PERIOD: ', period
    print p1
    print p2
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Days")
    plt.title("{}".format(name))

    #   plt.savefig('xycar_alldata.png')

    plt.show()

    return p2, order

#   want to do a weighted least squares fit

def fs_wtlsq_BIC(mjd, mag, name, errs, period=None, order=7):
    """
        This is the main fitting program I'm using at the moment for fourier series at least. It fits with least squares.
        Now try to add in an information criterion
    """
    #   initial guess for period
    ord = 1.
    if period==None:
        period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    p0 = [median, amplitude, 0.]
    
    #   least squares from python finds the parameters that minimize the error function, which is the difference between the model and the actual parameters
    p1, success = optimize.leastsq(weighted_error_series, p0, args=(mjd, mag, ord, period, errs))
    errors = np.sum(np.power(weighted_error_series(p1, mjd, mag, ord, period, errs), 2))
    add_order = True
    n = len(mjd)
    
    while (add_order == True) & (ord < order):
        ord2 = ord + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(weighted_error_series, p00, args=(mjd, mag, ord2, period, errs))
        errors2 = np.sum(np.power(weighted_error_series(p2, mjd, mag, ord2, period, errs), 2))
        #rederrors2 = errors2/
        k = 1 + ord*2.
        BIC = n*np.log(errors/n) + k*np.log(n)
        k2 = 1 + ord2*2.
        BIC2 = n*np.log(errors2/n) + k2*np.log(n)
        print "Errors", errors, errors2
        print "delta BIC", (BIC2-BIC)
        print "order", ord2
        if (BIC2-BIC) < 3:
            p1 = p00*1.
            errors = errors2*1.
            ord = ord2
        else:
            add_order = False
            ord = ord
    #   order = order - 1
    
    time = np.linspace(mjd.min(), mjd.max(), 300)
    #   time = np.linspace(53500, 53600, 300)
    plt.scatter(mjd, mag, alpha=0.4, color='grey')
    plt.plot(time, cosine_series(p2, time, ord, period), "r-")
    plt.gca().invert_yaxis()
    #   plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    #   plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[1],'\n' + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[0], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[0], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', ord
    print 'PERIOD: ', period
    print p1
    print p2
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Days")
    plt.title("{}".format(name))
    
    #   plt.savefig('xycar_alldata.png')
    
    plt.show()
    
    return p2, ord

def fs_wtlsq_BIC2(mjd, mag, name, errs, period=None, order=7):
    """
        This is the main fitting program I'm using at the moment for fourier series at least. It fits with least squares.
        Now try to add in an information criterion
        
        For some reason this version (which includes the part that extracts the best dates to fit with) doesn't work. So i'm going to just go with using the original BIC verson and writing a separate function for the extraction.
    """
    #   initial guess for period
    ord = 1.
    if period==None:
        period = get_period(name)

    ord = 1.
    period = get_period(cepheid)
    bin_width = period*4
    date_range = max(mjd)-min(mjd)
    nbins = np.around(date_range/bin_width)
    histout = np.histogram(mjd, nbins)
    binsizes = histout[0]
    print "binsizes", binsizes
    binedges = histout[1]
    bbin = np.where(binsizes == max(binsizes))[0]
    
    maxidx = np.argmin(np.abs(mjd - binedges[bbin+1]))
    minidx = np.argmin(np.abs(mjd - binedges[bbin]))
    print "bbin", bbin, "maxidx", maxidx, "minidx", minidx

    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    errs_tofit = errs[minidx:maxidx]
    mjd_tofit = mjd[minidx:maxidx]
    mag_tofit = mag[minidx:maxidx]
    
    p0 = [median, amplitude, 0.]
    
    #   least squares from python finds the parameters that minimize the error function, which is the difference between the model and the actual parameters
    p1, success = optimize.leastsq(weighted_error_series, p0, args=(mjd_tofit, mag_tofit, ord, period, errs_tofit))
    errors = np.sum(np.power(weighted_error_series(p1, mjd_tofit, mag_tofit, ord, period, errs_tofit), 2))
    add_order = True
    n = len(mjd_tofit)
    
    while (add_order == True) & (ord < order):
        ord2 = ord + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(weighted_error_series, p00, args=(mjd_tofit, mag_tofit, ord2, period, errs_tofit))
        errors2 = np.sum(np.power(weighted_error_series(p2, mjd_tofit, mag_tofit, ord2, period, errs_tofit), 2))
        #rederrors2 = errors2/
        k = 1 + ord*2.
        BIC = n*np.log(errors/n) + k*np.log(n)
        k2 = 1 + ord2*2.
        BIC2 = n*np.log(errors2/n) + k2*np.log(n)
        print "Errors", errors, errors2
        print "delta BIC", (BIC2-BIC)
        print "order", ord2
        if (BIC2-BIC) < 3:
            p1 = p00*1.
            errors = errors2*1.
            ord = ord2
        else:
            add_order = False
            ord = ord
    #   order = order - 1
    
    time = np.linspace(mjd_tofit.min(), mjd_tofit.max(), 300)
    #   time = np.linspace(53500, 53600, 300)
    plt.scatter(mjd_tofit, mag_tofit, alpha=0.4, color='grey')
    plt.plot(time, cosine_series(p2, time, ord, period), "r-")
    plt.gca().invert_yaxis()
    #   plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    #   plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[1],'\n' + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[0], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[0], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', ord
    print 'PERIOD: ', period
    print p1
    print p2
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Days")
    plt.title("{}".format(name))
    
    #   plt.savefig('xycar_alldata.png')
    
    plt.show()
    
    return p2, ord

def fs_wtlsq(mjd, mag, name, errs, period=None, order=7):
    #   initial guess for period
    ord = 1.
    if period==None:
        period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    p0 = [median, amplitude, 0.]
    
    #   least squares from python finds the parameters that minimize the error function, which is the difference between the model and the actual parameters
    p1, success = optimize.leastsq(weighted_error_series, p0, args=(mjd, mag, ord, period, errs))
    errors = np.sum(np.power(weighted_error_series(p1, mjd, mag, ord, period, errs), 2))
    add_order = True
    
    while (add_order == True) & (ord < order):
        ord = ord + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(weighted_error_series, p00, args=(mjd, mag, ord, period, errs))
        errors2 = np.sum(np.power(weighted_error_series(p2, mjd, mag, ord, period, errs), 2))
        if errors2 < errors:
            p1 = p00*1.
            errors = errors2*1.
        else:
            add_order = False
    #   order = order - 1
    
    time = np.linspace(mjd.min(), mjd.max(), 300)
    #   time = np.linspace(53500, 53600, 300)
    plt.scatter(mjd, mag, alpha=0.4, color='grey')
    plt.plot(time, cosine_series(p2, time, ord, period), "r-")
    plt.gca().invert_yaxis()
    #   plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    #   plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[1],'\n' + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[0], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[0], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', ord
    print 'PERIOD: ', period
    print p1
    print p2
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Days")
    plt.title("{}".format(name))
    
    #   plt.savefig('xycar_alldata.png')
    
    plt.show()
    
    return p2, order

def fs_wtlsq_lorentz(mjd, mag, name, errs, period=None, order=7):
    #   initial guess for period
    ord = 1.
    if period==None:
        period = get_period(name)
    median = approximate_median(mjd, mag)
    amplitude = approximate_amplitude(mjd, mag)
    zeroed_mag = mag - median
    
    p0 = [median, amplitude, 0.]
    
    #   least squares from python finds the parameters that minimize the error function, which is the difference between the model and the actual parameters
    p1, success = optimize.leastsq(weighted_error_series, p0, args=(mjd, mag, ord, period, errs))
    errors = np.sum(np.power(weighted_error_series(p1, mjd, mag, ord, period, errs), 2))
    add_order = True
    
    while (add_order == True) & (ord < order):
        ord = ord + 1
        p00 = np.concatenate([p1, [0.02, 0.]])
        p2, success = optimize.leastsq(weighted_error_series, p00, args=(mjd, mag, ord, period, errs))
        errors2 = np.sum(np.power(weighted_error_series(p2, mjd, mag, ord, period, errs), 2))
        if errors2 < errors:
            p1 = p00*1.
            errors = errors2*1.
        else:
            add_order = False
    #   order = order - 1
    
    time = np.linspace(mjd.min(), mjd.max(), 300)
    #   time = np.linspace(53500, 53600, 300)
    plt.scatter(mjd, mag, alpha=0.4, color='grey')
    plt.plot(time, cosine_series(p2, time, ord, period), "r-")
    plt.gca().invert_yaxis()
    #   plt.ylim([min(y), heapq.nlargest(2, y)[1]])
    #   plt.xlim([53500, 53600])
    print 'INITIAL GUESS' + '\n' + 'Amplitude:', p0[1],'\n' + 'Phase Shift: ', p0[2], '\n' + 'Offset: ', p0[0], '\n'
    print 'FIT PARAMETERS' + '\n' + 'Amplitude:', p1[1], '\n'  + 'Phase Shift: ', p1[2], '\n' + 'Offset: ', p1[0], '\n'
    print 'SUMMED ERRORS: ', errors
    print 'FIT ORDER: ', ord
    print 'PERIOD: ', period
    print p1
    print p2
    plt.ylabel("V-band Magnitude")
    plt.xlabel("Days")
    plt.title("{}".format(name))
    
    #   plt.savefig('xycar_alldata.png')
    
    plt.show()
    
    return p2, order

def best_fit_dates(mjd, period, nperiods=4., minpoints=100):
    """
        Picks out a small range of mjd over which to fit the cepheid. If there isn't a section with > 50 data points, it keeps making the bins bigger.
    
    """
    bin_width = period*nperiods
    date_range = max(mjd)-min(mjd)
    nbins = np.around(date_range/bin_width)
    histout = np.histogram(mjd, nbins)
    binsizes = histout[0]
    biggest_bin = max(binsizes)
    while(biggest_bin < minpoints):
        print "Not enough data points. Resizing bins."
        nperiods += 1
        bin_width = period*nperiods
        date_range = max(mjd)-min(mjd)
        nbins = np.around(date_range/bin_width)
        histout = np.histogram(mjd, nbins)
        binsizes = histout[0]
        biggest_bin = max(binsizes)
    
    binedges = histout[1]
    bbin = np.where(binsizes == max(binsizes))[0]
    print "binedges", binedges[bbin+1]
    print "binedges", binedges[bbin]
    
    maxidx = np.argmin(np.abs(mjd - binedges[bbin+1]))
    minidx = np.argmin(np.abs(mjd - binedges[bbin]))
    
    return bin_width, minidx, maxidx


"""
    class Fourier(object):
    
"""

if __name__=="__main__":
    
    import phase_folding as pf

    cepheid = 'aqpup'
    real_period = get_period(cepheid)

    ceph_data = get_sorted_data(cepheid)
    mjd = ceph_data[:,0]
    mag_v = ceph_data[:,1]
    epochs = ceph_data[:,2]
    phases = ceph_data[:,3]
    mag_v_err = ceph_data[:,4]

    good_ind = pf.simple_sigma_clip(mag_v, phases)
    good_ceph_dat = ceph_data[good_ind]
    good_mjd = mjd[good_ind]
    sortinds = np.argsort(good_mjd)
    sorted_good_ceph_dat = good_ceph_dat[sortinds]
    good_mjd = good_mjd[sortinds]
    good_mag_v = sorted_good_ceph_dat[:,1]
    good_epochs = sorted_good_ceph_dat[:,2]
    good_phases = sorted_good_ceph_dat[:,3]
    good_mag_v_err = sorted_good_ceph_dat[:,4]

    #coeff, order = fs_wtlsq(good_mjd, good_mag_v, cepheid, good_mag_v_err, order=12)

    binwidth, minidx, maxidx = best_fit_dates(good_mjd, period)
    mjd_tofit = mjd[minidx:maxidx]
    mag_tofit = mag_v[minidx:maxidx]
    phs_tofit = phases[minidx:maxidx]
    err_tofit = mag_v_err[minidx:maxidx]

    plt.clf()
    plt.scatter(phs_tofit, mag_tofit)
    plt.show()
    
    coeff2, order2 = fs_wtlsq_BIC(mjd_tofit, mag_tofit, cepheid, err_tofit, period=real_period, order=30)

    num_phase = 201 # this is just the number of phase points for the plot, won't really affect the run time
    phase_points = np.linspace(0, 1.0, num_phase)
    template_mags = cosine_series_phase(coeff, phase_points, order)
    template_mags2 = cosine_series_phase(coeff2, phase_points, order2)
    plt.clf
    plt.scatter(good_phases, good_mag_v, color='grey', alpha=0.3)
    #plt.plot(phase_points, template_mags, color='red', linewidth=2)
    plt.plot(phase_points, template_mags2, color='blue', linewidth=2)
    plt.ylabel("Magnitude")
    plt.xlabel("Phase")
    plt.show()


