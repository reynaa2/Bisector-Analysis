import numpy as np
from matplotlib import pyplot as plt
import math
from astropy.io import fits
import matplotlib.mlab as mlab

file = fits.open('apStarC-r5-2M16401588+3741030-55753.fits')
file.info()

print file[9]