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
    star_data = Cepheid_visual[number]
    name = star_data[0]

    IOMC_dat = get_IOMC_data(number)
    ASAS_dat = get_ASAS_data(name)

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
def get_ASAS_data(name):
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
        m2_asasiomc_old = asas_old['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
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

def get_IOMC_data(number):
    btime = 51544.5
    filename = "../iomc/IOMC_{}.fits".format(ids[number])
    hdulist = pyfits.open(filename)
    ceph_data = hdulist[1].data
    columns = hdulist[1].columns

    mag_v = ceph_data.field('mag_v')
    mjd = ceph_data.field('barytime') + btime
    errmag_v = ceph_data.field('errmag_v')
    exptime = ceph_data.field('exposure')
    pp = ceph_data.field('problems')

    period = per[number]
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

