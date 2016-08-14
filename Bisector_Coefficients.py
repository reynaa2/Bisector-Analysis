import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
from astropy.io import fits
import pylab as pl
import matplotlib.lines as mlines

# Primary issues: The file seems to stop writing out to the document at a particular 2MASSID
# Oddly the number of slope elements > DrewsList elements 

# Reading in the txt associated with the R and x-range values 
SB_x=np.loadtxt('SB9_Bisector_pts.csv', skiprows=1, dtype =float, usecols=[1,2,3])
DrewsList_x=np.loadtxt('DrewsList_Bisector_pts.lis', skiprows=1, dtype=float, usecols=[1,2,3])
DrewsList2_x=np.loadtxt('DrewsList2_Bisector_pts.lis', skiprows=1, dtype=float, usecols=[1,2,3])
DR12_random_x=np.loadtxt('APOGEE_Random_Bisector_pts.lis',skiprows=1, dtype=float, usecols=[1,2,3] )

MASSID=np.loadtxt('DrewsList_Bisector_pts.lis',dtype = 'str', usecols=[0])
twoMASSID = MASSID[1:]
print(len(twoMASSID))
y = [0.1,0.3,0.5,0.7,0.9]
y1 = [0.1,0.3,0.5]
y2 = [0.1,0.3,0.5,0.7]

#Plotting bisectors of only three points
plt.figure(figsize=(8,8))
slope=[]
poly_coeff=[]
for i in range(len(DrewsList_x)):
    
    plt.figure(figsize=(8,8))
    plt.plot(DrewsList_x[i],y1,'ro')
    plt.xlabel('CCF lag')
    plt.ylabel('Agreement')
    plt.title('Bisectors for '+str(twoMASSID[i]))
    plt.plot(DrewsList_x[i], np.poly1d(np.polyfit(DrewsList_x[i], y1, 1))(DrewsList_x[i]));
    slope.append(np.polyfit(DrewsList_x[i],y1,1));
    m = np.polyfit(DrewsList_x[i],y1,1);

    plt.plot(DrewsList_x[i],np.poly1d(np.polyfit(DrewsList_x[i],y1,4))(DrewsList_x[i]));
    poly_coeff.append(np.poly1d(np.polyfit(DrewsList_x[i],y1,4)));
    #plt.savefig('Bisector for '+str(twoMASSID[i])+' with 3pts')
    #plt.show()
    
#Writing a csv to hold the 2MASSID and bisector linear coefficients for three points
outfile2 = open('DrewsList_Coefficients_for_3pts.csv','w')
for i in range(len(DrewsList_x)):
    outfile2.write(str(twoMASSID[i])+'\t'+str(slope[i][0])+'\n')

outfile2.close()
