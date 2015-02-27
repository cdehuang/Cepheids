#
#   chuang
#
#   02.20.15
#
#

"""
    Code to figure out what phase the Cepheid's light curve will be at a given time (in Julian Days). Input parameters should be the light curve information and the Julian Date at the date of observation.
    Incorporates IDL code for the Yoachim et al. templates that I got from Samantha so it requires pIDLy
"""

#   python packages to import
import numpy as np
import pidly
import astropy.io.ascii as ascii
import asciitable

#   python modules to import (that were written by me)
import ft
import pdm

def findphs(lc, mjd):
    """
        lc: location of the Cepheid lightcurve file
        mjd: date of observation for which we want the phase.
        
        Going to try to avoid having multiple return statements so this program will basically just call all the other ones. So the first thing this one will have to do is determine if there is a period change for the Cepheid (if there isn't, then it's going to be really easy to get a period.
        
        Not sure yet if it is a good idea to write this as a function or if it should just be a python script
    """
    return "TO DO"

#initialize IDL enviroment in python
idl = pidly.IDL()

#light curve location/make light curve
cepheid_name = "xycar"
fname = 'data'
file = ft.make_data_file(cepheid_name, fname)
lc = asciitable.read(file)


#get the period
period = ft.get_period(cepheid_name)

#make a template Cepheid with the appropriate period and see how it works. Maybe here try phase-dispersion analysis, and if the phase is too disperse, just
per, theta = pdm.find_period(file, period)

#either way, we first have to fit the Cepheid light curve with a template.
#We will use the Yoachim et al templates and Cepheid-fitting code obtained from online
mjd = lc['mjd']
fitparams = []
idl('print', 'hello')

#fit the Cepheid light curve for that period to the data and see how good of a job it does (it will fit the starting phase, the amplitude, and the period)
#maybe here it is best to start off with the Fourier series fitting (to get quick estimates) and then try to fit an actual light curve to it. Or maybe it would be faster just to fit light curves to it, I don't really know.


#what the tolerance is for not fitting. I think something like 0.1 is pretty realistic. So if this is greater, then assume that the period is actually changing and try out various period changes on it. This is where we actually care about fitting the Cepheid light curve to the
threshold = 0.1
if (theta >= threshold):
    print "The Cepheid's period may be changing."


if __name__=="__main__":
#initialize IDL environment
    idl = pidly.IDL()

