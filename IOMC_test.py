#   chuang
#   12.8.2014

import numpy as np
import astropy.io.ascii as ascii
import pyfits
import matplotlib.pyplot as plt
import asciitable
from scipy.interpolate import UnivariateSpline
import os.path
from matplotlib import cm
from numpy import *

#
#
#   Change so that you can just type in the name of the Cepheid instead of the number.
#   Also change the whole thing so that it's a function instead. Then you can go and go and make it so that if you don't give it an input period, you can guess one.
Cepheid_visual = np.genfromtxt('../iomc/Cepheid_visual.txt', dtype=None, names=['name', 'id', 'cr1', 'cr2', 'ra1', 'ra2', 'ra3', 'dec1', 'dec2', 'dec3', 'radec', 'decdec', 'cr11', 'type', 'cr22', 'cr33', 'm1', 'err', 'mmax', 'mmin', 'per', 'gal_l', 'gal_b', 'Schafly Av', 'SFD Av', 'mag/kpc'], skip_header=2, filling_values=-1.0)
ids = Cepheid_visual['id']
names = Cepheid_visual['name']
id_numbers = np.arange(len(ids))
cepheids = {}
for i in id_numbers:
    cepheids[names[i]] = i
per = Cepheid_visual['per']

name = 'svul'
number = cepheids[name]
star_data = Cepheid_visual[number]

#
#   INTEGRAL-OMC data
#
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

#Getting the smallest magnitude for plotting purposes and throwing out the worst point
minval = np.min(mag_v)
mi = np.min(n for n in mag_v if n!=minval)
maxval = np.max(mag_v)
ma = np.max(n for n in mag_v if n!=maxval)

period = per[number]
newjd = mjd - np.min(mjd)

phase = mjd/period - np.fix(mjd/period)
sortidxs = np.argsort(phase)
sorted_phase = phase[sortidxs]
sorted_mag_v = mag_v[sortidxs]
sorted_mjd = mjd[sortidxs]
arr = np.arange(len(ceph_data))
gd = arr[(exptime[arr] < 200) & (pp[arr] == np.median(pp))]
mag_v = mag_v[gd]
mjd = mjd[gd]
phase = phase[gd]
epochs = np.fix(mjd/period)

t = mjd
y = mag_v
ep = epochs
phs = phase

#   Try to bin data that is less separated than an hour into each epoch and extract the one median point from each binning. < 0.04 days

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

phn = binned_t/period - np.fix(binned_t/period)
sort_idxs = np.argsort(phn)
sort_phase = phn[sort_idxs]
sort_mag = binned_mag[sort_idxs]

double_phase = np.concatenate([phs, phs+1.0], 1)
double_mag = np.concatenate([mag_v, mag_v], 1)

"""

plt_phase_old = old_phases - phaseshift
    opphase = []
    for i in plt_phase_old:
        if i < 0.0:
            opphase.append(i + 1.0)
        else:
            opphase.append(i)
    opphase2 = [i + 1.0 for i in opphase]
    double_old_phase = [opphase, opphase2]
    double_old_mag = [old_m2, old_m2]


"""
n_bins = (np.max(t) - np.min(t))*24
count, bin = np.histogram(t, n_bins)
bins = np.linspace(np.min(t), np.max(t), n_bins)
bin_ind = np.digitize(t, bins)

digitize_t = []
digitize_mag = []

segs_end_ind = np.concatenate([np.where(np.diff(bin_ind) != 0)[0] + 1, [len(bin_ind)]], 1)
segs_begin_ind = np.concatenate([[0], np.where(np.diff(bin_ind) != 0 )[0] + 1],1)
digitize_num = np.arange(len(segs_end_ind))
for i in digitize_num:
    time_segment = t[segs_begin_ind[i]:segs_end_ind[i]]
    mag_segment = y[segs_begin_ind[i]:segs_end_ind[i]]
    phase_segment = time_segment/period - np.fix(time_segment/period)
    print mag_segment
    median_time = np.median(time_segment)
    median_mag = np.median(mag_segment)
    digitize_t.append(median_time)
    digitize_mag.append(median_mag)

digitize_t = np.asarray(digitize_t)
digitize_mag = np.asarray(digitize_mag)
digiph = digitize_t/period - np.fix(digitize_t/period)
#for i in np.unique(bin_ind):

#
#   Old ASAS data (the set that you can download from their website)
#
old_asas_file = "../asas_old/{}.txt".format(star_data[0])
if os.path.isfile(old_asas_file) != False:
    asas_old = np.genfromtxt(old_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
    asasdate_old = asas_old['asasdate'] + 50000.0
    old_asasphasenum = np.fix(asasdate_old/period)
    old_asasphase = asasdate_old/period - np.fix(asasdate_old/period)
    grade_old = asas_old['grade']
    old_arr2 = np.arange(len(asas_old))
    old_asasgd = old_arr2[grade_old[old_arr2] == 'A']
    m2_old = asas_old['m2']
    m2_asasiomc_old = asas_old['m2'] - 0.03 #to make the data agree with the v-band magnitude from IOMC
    old_phases = old_asasphase[old_asasgd]
    old_m2 = m2_old[old_asasgd]
    s = UnivariateSpline(old_phases, old_m2)
    all_phasenumbers = np.concatenate([old_asasphasenum[old_asasgd], ep],1)
    all_magnitudes = np.concatenate([old_m2, y], 1)
    all_dates = np.concatenate([asasdate_old[old_asasgd], t], 1)
    all_phases = np.concatenate([old_phases, phs], 1)
    deltaphases = np.linspace(0, 1, 101)
    phasemags = s(deltaphases)
    phaseshift = deltaphases[np.argmax(phasemags)]
    plt_phase = sorted_phase - phaseshift
    pphase = []
    for i in plt_phase:
        if i < 0.0:
            pphase.append(i + 1.0)
        else:
            pphase.append(i)
    pphase2 = [i + 1.0 for i in pphase]
    double_phase = [pphase, pphase2]
    double_mag = [sorted_mag_v, sorted_mag_v]

    plt_phase_old = old_phases - phaseshift
    opphase = []
    for i in plt_phase_old:
        if i < 0.0:
            opphase.append(i + 1.0)
        else:
            opphase.append(i)
    opphase2 = [i + 1.0 for i in opphase]
    double_old_phase = [opphase, opphase2]
    double_old_mag = [old_m2, old_m2]
    print "old asas -- yes"

#
#   New ASAS data (the set that Adam got from Dorota by email correspondence)
#
cepheid_filename = '_'.join([name[:-3],name[-3:]])
new_asas_file = "../new_asas/{}.txt".format(cepheid_filename)
if os.path.isfile(new_asas_file) != False:
    asas_new = np.genfromtxt(new_asas_file, dtype=None, names = ['asasdate', 'm2', 'm0', 'm1', 'm3', 'm4', 'e2', 'e0', 'e1', 'e3', 'e4', 'grade', 'frame'])
    new_asasdate = asas_new['asasdate'] + 50000.0
    #    new_asasdate2 = new_asasdate + 0.23
    new_asasphasenum = np.fix(new_asasdate/period)
    new_asasphase = new_asasdate/period - np.fix(new_asasdate/period)
    #    new_asasphase2 = new_asasdate2/period - np.fix(new_asasdate2/period)
    grade_new = asas_new['grade']
    new_arr2 =  np.arange(len(asas_new))
    new_asasgd = new_arr2[grade_new[new_arr2] == 'A']
    m2_new = asas_new['m2']
    m2_asasiomc_new = asas_new['m2'] - 0.03
    new_m2 = m2_new[new_asasgd] + 0.01
    new_phases = new_asasphase[new_asasgd]

    plt_phase_new = new_phases - phaseshift
    npphase = []
    for i in plt_phase_new:
        if i < 0.0:
            npphase.append(i + 1.0)
        else:
            npphase.append(i)
    npphase2 = [i + 1.0 for i in npphase]
    double_new_phase = [npphase, npphase2]
    double_new_mag = [new_m2, new_m2]
    all_phasenumbers = np.concatenate([new_asasphasenum[new_asasgd], all_phasenumbers],1)
    all_magnitudes = np.concatenate([new_m2,all_magnitudes], 1)
    all_dates = np.concatenate([new_asasdate[new_asasgd], all_dates],1)
    all_phases = np.concatenate([new_phases,all_phases], 1)
    print "new_asas -- yes"

plt.clf()
plt.title("{}".format(star_data[0]))
plt.xlabel("Phase")
plt.ylabel("V-Band Magnitude")
if os.path.isfile(old_asas_file) != False:
    dap = np.concatenate([all_phases, all_phases+1.0], 1)
    dam = np.concatenate([all_magnitudes, all_magnitudes], 1)
    plt.scatter(dap, dam, alpha=0.6, color='grey')
else:
    plt.scatter(double_phase, double_mag, alpha=0.5, color='grey')
plt.scatter(sort_phase, sort_mag, alpha=0.6, color='red')
plt.ylim([minval-0.15, maxval+0.15])
plt.gca().invert_yaxis()
plt.show()

"""

svul_dates = [56626.26084,56745.26084,56957.26084]
med_mag = np.median(all_magnitudes)
svul_datesm = [med_mag, med_mag, med_mag]

plt.clf()
plt.scatter(svul_dates, svul_datesm, alpha=0.85, color='blue')
plt.scatter(all_dates, all_magnitudes, color='grey', alpha=0.7)
plt.title("S Vul Observations")
plt.xlabel("V-band Magnitude")
plt.ylabel("Julian Date")
plt.ylim([minval-0.15, maxval+0.15])
plt.gca().invert_yaxis()
plt.savefig("../output/svul_obs.png")
plt.show()
"""
