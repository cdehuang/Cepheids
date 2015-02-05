#   chuang
#   12.4.2014

import numpy as np
import astropy.io.ascii as ascii
import pyfits
import matplotlib.pyplot as plt
import asciitable
from scipy.interpolate import UnivariateSpline
import os.path
from matplotlib import cm
from scipy import *
#import scipy.optimize as op
from scipy import optimize
import george
from george import kernels
import emcee

#   Take the data, try to do a discrete fourier transform on it. The Integral data seems to not be well-behaved. Try it without integral

#
#   Text file listing all of the Cepheids and their periods, etc.
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

name = 'xzcar'
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
minval = min(mag_v)
mi = min(n for n in mag_v if n!=minval)
maxval = max(mag_v)
ma = max(n for n in mag_v if n!=maxval)

period = 16.653
newjd = mjd - min(mjd)

phase = mjd/period - np.fix(mjd/period)
sortidxs = np.argsort(phase)
sorted_phase = phase[sortidxs]
sorted_mag_v = mag_v[sortidxs]
arr = np.arange(len(ceph_data))
gd = arr[(exptime[arr] < 200) & (pp[arr] == np.median(pp))]
phases = np.fix(mjd[gd]/period)

all_phasenumbers = phases*1.0
all_magnitudes = mag_v*1.0
all_dates = mjd[gd]*1.0
all_phases = phase*1.0

asasonly_dates = []
asasonly_magnitudes = []

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
    all_phasenumbers = np.concatenate([old_asasphasenum[old_asasgd], all_phasenumbers],1)
    all_magnitudes = np.concatenate([old_m2, all_magnitudes], 1)
    all_dates = np.concatenate([asasdate_old[old_asasgd], all_dates], 1)
    all_phases = np.concatenate([old_phases, all_phases], 1)
    asasonly_dates = np.concatenate([asasdate_old[old_asasgd], asasonly_dates], 1)
    asasonly_magnitudes = np.concatenate([old_m2, asasonly_magnitudes], 1)
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
    asasonly_dates = np.concatenate([new_asasdate[new_asasgd], asasonly_dates], 1)
    asasonly_magnitudes = np.concatenate([new_m2, asasonly_magnitudes], 1)
    print "new_asas -- yes"

#all_phasenumbers = np.concatenate([old_asasphasenum[old_asasgd], new_asasphasenum[new_asasgd], phases],1)
#all_magnitudes = np.concatenate([old_m2, new_m2, mag_v], 1)
#all_dates = np.concatenate([asasdate_old[old_asasgd], new_asasdate[new_asasgd], mjd[gd]], 1)
#all_phases = np.concatenate([old_phases, new_phases, phase], 1)

####
#   Attempts at fitting things
####


class Parameter:
    def __init__(self, value):
        self.value = value
    
    def set(self, value):
        self.value = value
    
    def __call__(self):
        return self.value

def fit(function, parameters, y, x = None):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)
    
    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    optimize.leastsq(f, p)

#   now give the initial parameters

#   sort the data points, mind the spacing between them as well. Try picking out some data points with separation greater than some number and all other data points with separation smaller than some other number. or maybe throw out the IOMC data when it comes to fitting?
#   Alternatively, I could try writing a program that folds the cepheid light curves and then calculates the dispersion (in time) and tries to minimize it.
sortinds = np.argsort(all_dates)
t = all_dates[sortinds]
y = all_magnitudes[sortinds]
sortinds2 = np.argsort(asasonly_dates)
t2 = asasonly_dates[sortinds2]
y2 = asasonly_magnitudes[sortinds2]

#   define the function:
fitfunc = lambda p, x: p[0]*cos(2*pi/p[1]*x+p[2]) + p[3]#Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y
p0 = [0.4, 16.6, 0, 8.5]
p1, success = optimize.leastsq(errfunc, p0[:], args=(t, y))
p2, success = optimize.leastsq(errfunc, p0[:], args=(t2, y2))

print p1

t_seg = t2[200:300]
y_seg = y2[200:300]
time = np.linspace(min(t_seg), max(t_seg), 200)

plt.clf()
plt.scatter(t_seg, y_seg)
plt.plot(time, fitfunc(p2, time))
plt.plot(time, fitfunc(p0, time), color='green')
plt.show()


#   Find largest "consecutive run
run_ends = np.where(np.diff(t) > 10)[0] + 1
max_run_end_loc = np.where(np.diff(run_ends)>10)[0].max() + 1
max_run_end_ind = run_ends[max_run_end_loc]
max_run_begin_loc = max_run_end_loc - 1
max_run_begin_ind = run_ends[max_run_begin_loc]
t_segment = t[max_run_begin_ind: max_run_end_ind]
y_segment = y[max_run_begin_ind: max_run_end_ind]

kernel = kernels.CosineKernel(16.6) + kernels.WhiteKernel(2)
gp = george.GP(kernel)
gp.compute(t_seg)
print(gp.lnlikelihood(y_seg))
print(gp.grad_lnlikelihood(y_seg))

def nll(p):
    gp.kernel[:] = p
    ll = gp.lnlikelihood(y, quiet=True)

    return -ll if np.isinfinite(ll) else 1E25

def grad_nll(p):
    gp.kernel[:] = p
    return -gp.grad_lnlikelihood(y, quiet=True)

gp.compute(t_seg)

print(gp.lnlikelihood(y_seg))
p0 = gp.kernel.vector
results = op.minimize(nll, p0, jac=grad_nll)

gp.kernel[:] = results.x
print(gp.lnlikelihood(y_seg))

"""
def lnprob(p):
    # Trivial improper prior: uniform in the log.
    if np.any((-10 > p) + (p > 10)):
        return - np.inf
    lnprior = 0.0

    # Update the kernel and compute the likelihood.
    kernel.pars = np.exp(p)
    return lnprior + gp.likelihood(y, quiet=True)

# You need to compute the GP once before starting. Then the sample list will be saved.
gp.compute(t)

#Set up sampler
nwalkers, ndim = 36, len(kernel)
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)

#Initialize walkers.
p0 = [np.log(kernel.pars) + 1E-4*np.random.randn(ndim) for i in range(nwalkers)]

print("Running burn-in")
p0, _, _ = sampler.run_mcmc(p0, 2000)

print("Running production chain")
sampler.run_mcmc(p0, 2000)
"""
"""
####
#   Plotting everything
####
plt.clf()
plt.title("{}".format(star_data[0]))
plt.xlabel("Phase")
plt.ylabel("V-Band Magnitude")
"""
"""
if os.path.isfile(new_asas_file) != False:
    plt.scatter(double_new_phase, double_new_mag, color='teal', alpha=0.7, edgecolors='none')
"""
"""
period_list = set(all_phasenumbers)
x = 0
for i in period_list:
    period_date_idxs = np.asarray(np.where(all_phasenumbers == i))
    period_mags = all_magnitudes[period_date_idxs]
    period_dates = all_dates[period_date_idxs]
    period_phases = all_phases[period_date_idxs]
    period_phases2 = [i + 1.0 for i in period_phases]
    double_period_phase = [period_phases, period_phases2]
    double_period_mag = [period_mags, period_mags]
    plt.scatter(double_period_phase, double_period_mag, color=cm.jet(1.*x/len(period_list)), alpha=0.7, edgecolors='none')
    x +=1
#for i in len(all_phasenumbers):
#plt.scatter(double_phase, double_mag, color='pink', alpha=0.85, edgecolors='none')
plt.xlim([-0.2, 2.2])
if os.path.isfile(old_asas_file) != False:
    plt.ylim([max(old_m2) +0.3 , min(old_m2)-0.3])
    plt.text(-0.1, min(old_m2)-0.15, "period = {} days".format(period), fontsize=14)
elif os.path.isfile(new_asas_file) != False:
    plt.ylim([max(new_m2)+0.3, min(new_m2)-0.3])
    plt.text(-0.1, min(new_m2)-0.15, "period = {} days".format(period), fontsize=14)
else:
    plt.ylim([max(mag_v) + 0.3, min(mag_v) - 0.3])
    plt.text(-0.1, mi-0.15, "period = {} days".format(period), fontsize=14)
plt.savefig("../output/{}_periodfs.png".format(name))
plt.show()
"""
