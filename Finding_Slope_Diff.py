import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
from astropy.io import fits
import pylab as pl
import matplotlib.lines as mlines

SB_x=np.loadtxt('SB9_Coefficients_for_3pts.csv',dtype =float, usecols=[1])
twoMASSID = np.loadtxt('SB9_Coefficients_for_3pts.csv',dtype =str, usecols=[0])

DrewsList_x=np.loadtxt('DrewsList_Coefficients_for_3pts.csv', skiprows=1, dtype=float, usecols=[1])
TwoMASSID = np.loadtxt('DrewsList_Coefficients_for_3pts.csv', skiprows=1, dtype=str, usecols=[0])

IDs = []
for j in range(SB_x.size):
    if str(twoMASSID[j:]) == str(twoMASSID[j]):
         IDs[j] = twoMASSID[j]
         print IDs
#DrewsList2_x=np.loadtxt('DrewsList2_Coefficients_for_3pts.csv', skiprows=1, dtype=float, usecols=[1])
#DR12_random_x=np.loadtxt('APOGEE_random_Coefficients_for_3pts.csv',skiprows=1, dtype=float, usecols=[1] )


#Finding the difference in slope per each visit
diff_in_x=[]
for i in range(SB_x.size):
    for i in range(len(SB_x)):
        diff_in_x = diff_in_x + [SB_x[i:]-j for j in SB_x[i+1:]] 
        print(diff_in_x)
    
print(len(diff_in_x))
            #difm=difm+[ccm[i]-j for j in ccm[i+1:]]
#print(len(diff_in_x))