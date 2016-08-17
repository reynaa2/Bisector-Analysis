import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
from astropy.io import fits
import pylab as pl
import matplotlib.lines as mlines
from scipy.interpolate import interp1d


def measure_r(comb,velocity_axis,peak_loc):
    n_cc = comb.size
    ccf_max = max(comb,peak_loc)

    #Split the cross correlation array on its peak 
    if n_cc > peak_loc*2:
        low_comb = comb[0:peak_loc]
        high_comb = comb[peak_loc: 2*peak_loc]
    else:
            low_comb = comb[2*peak_loc-n_cc+1:peak_loc]
            high_comb = comb[peak_loc:n_cc-1]
        
    #Flip the low side to get asymmetric component
    low_comb = np.fliplr(low_comb)

    if low_comb.size == high_comb.size:
        one_over_2n = 1./(2.*low_comb.size)
        sigma_a = math.sqrt(one_over_2n* np.sum((high_comb - low_comb)^2.0)
        height  = comb[peak_loc]
        r = height/(math.sqrt(2.0)*sigma_a)
        print(r)
return r

def ccf_width(ccf, ccflag):
    max_lag = max(ccflag)
    min_lag = min(ccflag)
    ccf_max = max(ccf)
    ccf_min = min(ccf)
    height = ccf_max - ccf_min
    above_half = while(ccf >= 0.5*ccf_max, n_above_half)
    if  n_above_half > 1.0:
         low_above_half = min(ccflag[above_half])
         high_above_half = max(ccflag[above_half])
         fwhm = high_above_half - low_above_half
    elif:
        fwhw = 1.0
return fwhm

def CCF_height(ccfm,ccflag):
    max_lag = max(ccflag)
    min_lag = min(ccflag)
    ccf_max = max(ccf)
    ccf_min = min(ccf)
    ccf_height = ccf_max - ccf_min
return ccf_height

def threshold(ccf,ccflag):
     max_lag = max(ccflag)
     min_lag = min(ccflag)
     ccf_max = max(ccf)
     ccf_min = min(ccf)
     height = ccf_max - ccf_min
    #Normalize the CCF in the y-direction so peak=1 and floor=0
     norm_ccf = (ccf - ccf_min)/(ccf_max-ccf_min)

   #Make finer grid via interpolation. Create grid size
    finer_grid = np.arange((max_lag-min_lag)*20.)/20. + min_lag
    finer_ccf = interpld(norm_ccf, ccflag,finer_grid)
  
  #Find the ccf threshold
   ccf_median = median(finer_ccf)
   threshold = 2.0 * ccf_median
return threshold

def measure_bisector(ccflag, ccf):
    #Find height of the CCF 
    ccf_max = max(ccf)
    ccf_min = min(ccf)
    height = ccf_max-ccf_min
    #Find boundaries that split the CCF into n-zones
    n_boundaires=6
    zone_boundaries = ccf_min+(height*findgen(n_boundaires)/(n_boundaires-1))
    #Make arrays to store teh vertical middle of each zone and the resultant horizontal bisector [note that we can only measure n-zones -1 bisectors, as there are more boundaries]
    bisector_postions =  np.zeros((2, n_boundaries-1)
    #Loop through the zones and find bisector pt of each
    for i in range(n_boundaries-1):
        #save y-position of the bisector
        bisector_positions[1,i]=(zone_boundaries[i]+zone_boundaires[1+i])/2.0
        #Identify portion of CCF inbetween boundaries
        in_boundaries = while(ccf >= zone_boundaries[i] and ccf <= zone_boundaries[i+1], n_in_boundaries)

        #Find the mean x position between these boundaries, it it exists. If not save -9999 to indicate measurement not real
        if n_boundaries >= 1:
            bisector_positions[0,i] = mean(ccflag[in_boundaries])
        elif:
            bisector_positions[0,i] = -9999
return bisector_positions

def measure_apogee_CCF_bisector(ccf, ccflag):
    max_lag=max(ccflag)
    min_lag=min(ccflag)
    ccf_max=max(ccf)
    ccf_min=min(ccf)
    height = ccf_max-ccf_min
    normalized_ccf = (ccf-ccf_min)/(ccf_max-ccf_min)
    #Make a grid with 10x sampling to be less suseptible to sampling effects
    finer_grid = np.arange((max_lag-min_lag)*20.)/20.+min_lag
    #Interpolate the CCF into this finer grid
    finer_ccf = interp1d(normalized_ccf,ccflag,finer_grid)
    #Find CCF's base level
    ccf_median = median(finer_ccf)
    truncate = 0.0
    to_analyze = while(finer_ccf > tuncate, n_analyze)
    if n_analyze > 3:
        bisector = measure_bisector(finer_grid[to_analyze-2], finer_ccf[to_analyze-2])
    #zoom in to the part of the CCF we are analyzing; find the plot
    #range we need to see that clearly.
    min_lag_to_analyze = min(finer_grid[to_analyze])
    max_lag_to_analyze = max(finer_grid[to_analyze])
    zoom_range = max_lag_to_analyze - min_lag_to_analyze
    #LOADCT, 0, /SILENT
    #PLOT, finer_grid, finer_ccf, XRANGE = [min_lag_to_analyze-0.1*zoom_range, max_lag_to_analyze+0.1*zoom_range], /XSTY , PSYM=10  ;, COLOR = 255
    #LOADCT, 13, /SILENT
    #OPLOT, [-250,250], [ccf_median, ccf_median], COLOR = 255
    #OPLOT, [-250,250], [ccf_median, ccf_median], COLOR = 80

return bisector

# Read in list holding apStar and apStarC files
twoMASSID, field, path = np.loadtxt('DR12_SB9_SB2s.lis',skiprows=1,usecols=[4,5,7])
TwoMassID, Field, Path = np.loadtxt('DrewsList_BC1.lis',skiprows=1,usecols=[0,1,2])
massID, apogee_loc, path_loc = np.loadtxt('DrewsList_BC2.lis',skiprows=1,usecols=[0,1,2])















