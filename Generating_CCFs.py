import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
import os
from astropy.io import fits

DR12 = 'C:\Users\Adela\Documents\Desktop\Research_2016\APOGEE_DR12_ApStar'
DR12_subdirectories = os.listdir(DR12)

ccf = fits.open('apStarC-r5-2M21432388+4209512.fits')
#print ccf.info()
ccm = ccf[9].data['CCF'][0][0]
spectra = ccf[1].data
#print spectra
ccf1= ccf[9].data['CCF'][0].shape
#print ccf1


