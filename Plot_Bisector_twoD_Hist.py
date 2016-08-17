import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
from astropy.io import fits
import pylab as pl
import matplotlib.lines as mlines

# Function developed to find the index of values in the most populated bin of 2D histogram (Developed by Graham)
def find_index(x_val, y_val, x_range, y_range, x_ind, y_ind):
    xmin=x_range[x_ind]
    xmax=x_range[x_ind+1]
    ymin=y_range[y_ind]
    ymax=y_range[y_ind+1]
    x_true=np.where(np.logical_and(x_val>=xmin,x_val<xmax))[0]
    y_true=np.where(np.logical_and(y_val>=ymin, y_val<ymax))[0]
    true_indicies=[i for i in x_true if i in y_true]
    return true_indicies #Returns the indicies for the 2MASSID 

# Reading in the lists associated with the R and x-range values 
SBR, SBXmax = np.loadtxt('SB9_Rmin_X_Max.lis',usecols=[2,3], skiprows=1,unpack=True, dtype=float)

# Drew's List Part 1 (contains 50 candidates)
DLR, DLXmax = np.loadtxt('Drews_List_Rmin_XR_Max.lis',usecols=[2,3], skiprows=1,unpack=True, dtype=float)

# Drew's List Part 2 (contains 1158 candidates)
DLR2, DLXmax2 = np.loadtxt('Drews_List2_Rmin_XR_Max.lis',usecols=[2,3], skiprows=1,unpack=True, dtype=float)

# Controlled sample of 713 random APOGEE-DR12 stars that contain non-binaries
ALR, ALX = np.loadtxt('APOGEE_Rand_Rmin_XR_Max.lis',usecols=[2,3], skiprows=1,unpack=True, dtype=float)

# Mass IDs for all lists
massID = np.loadtxt('SB9_Rmin_X_Max.lis',skiprows=1,dtype = str, usecols=[1])
MassID = np.loadtxt('Drews_List_Rmin_XR_Max.lis',skiprows=1,dtype=str,usecols=[1])
DmassID = np.loadtxt('Drews_List2_Rmin_XR_Max.lis',skiprows=1,dtype=str,usecols=[1])
AmassID = np.loadtxt('APOGEE_Rand_Rmin_XR_Max.lis',skiprows=1,dtype=str,usecols=[1])
# Field IDs for all lists
SBlocationID = np.loadtxt('SB9_Rmin_X_Max.lis',skiprows=1,dtype = str, usecols=[0])
DL_locationID = np.loadtxt('Drews_List_Rmin_XR_Max.lis',skiprows=1,dtype = str, usecols=[0])
DL2_locationID = np.loadtxt('Drews_List2_Rmin_XR_Max.lis',skiprows=1,dtype=str,usecols=[0])
APlocationID = np.loadtxt('APOGEE_Rand_Rmin_XR_Max.lis',skiprows=1,dtype=str,usecols=[0])

# All the mass IDs from the binary lists.
IDs = np.append(massID, np.append(MassID, DmassID))
# All mass IDs for the non-binary and binary list
comb=np.append(massID,np.append(MassID, DmassID))
MassIDs = np.append(comb,AmassID)


# All combined location IDs
field = np.append(SBlocationID, np.append(DL_locationID, DL2_locationID))
locationID = np.append(field, APlocationID)

#For R min vs max x-range of binary and controlled samples
outfile1=open('Combined Candidate Log Values.csv', 'w')
xs = np.append(SBXmax ,np.append(DLXmax, DLXmax2))
ys = np.append( SBR ,np.append(DLR, DLR2))
BinX = np.append(xs, ALX)
BinR = np.append(ys ,ALR)

xmax=[]
R=[]
for j in range(len(BinX)):
    xmax.append(math.log10(BinX[j]))
    R.append(math.log10(BinR[j]))
    outfile1.write(str(locationID[j])+'\t'+ str(MassIDs[j])+'\t'+str(xmax[j])+'\t'+str(R[j])+ '\t' +str(BinX[j])+ '\t'+ str(BinR[j])+'\n')
outfile1.close()

perBin = 0.025
binCt = int(perBin*len(xmax))
Grid1, x_range, y_range =np.histogram2d(xmax,R, bins=binCt, range=None, normed=False, weights=None)
arg  = np.argmax(Grid1)
t1 = arg%Grid1.shape[0]
t2 = arg/Grid1.shape[0]
found_ind=find_index(xmax, R, x_range, y_range, t2, t1)
print("Graham's algorithm found {0} indicies".format(found_ind))

#Call the index of the function's generated value for found_ind
print(MassIDs[found_ind])

plt.figure(figsize=(8,6))
CS = plt.pcolor(x_range, y_range,Grid1, cmap='RdBu')
plt.colorbar(label = 'Number of values')
plt.xlabel('log(max x-ranges)')
plt.ylabel('log(R min)')
plt.title('Combined Candidates Binary and Non-Binary' )
plt.savefig('Combined Cand. 2DHist for Log Values.jpg', bbox_inches='tight') 
#--------------------------------------------------
#For R min vs max x-range (for binary sample only-2011)
outfile=open('Binary Log Values.csv', 'w')
x = np.append(SBXmax, np.append(DLXmax, DLXmax2))
y = np.append(SBR, np.append(DLR,DLR2))

Xmax=[]
Rmin=[]
for i in range(len(x)):
    Xmax.append(math.log10(x[i]))
    Rmin.append(math.log10(y[i]))
    outfile.write( str(field[i])+'\t'+str(IDs[i])+'\t'+str(Xmax[i])+'\t'+str(Rmin[i])+ '\t' +str(x[i])+ '\t'+ str(y[i])+'\n')
outfile.close()

perBin = 0.04
binCt = int(perBin*len(Xmax))
Grid, X1_range, Y1_range =np.histogram2d(Xmax, Rmin, bins=binCt, range=None, normed=False, weights=None)

plt.figure(figsize=(8,6))
C = plt.pcolor(X1_range, Y1_range,Grid, cmap='RdBu')
plt.colorbar(label = 'Number of values')
plt.xlabel('log(max x-ranges)')
plt.ylabel('log(R min)')
plt.title('Combined Candidates Binary')
plt.savefig('Combined Binary Cand. Log Values.jpg', bbox_inches='tight')
max_grid  = np.argmax(Grid)
ta = max_grid%Grid.shape[0]
tb = max_grid/Grid.shape[0]
found_index=find_index(Xmax, Rmin, X1_range, Y1_range, tb, ta)
print("Graham's algorithm found {0} indicies".format(found_index))
print(IDs[found_index])

# Find the largest and smallest x-range max value 
XRmax = np.loadtxt('Binary Log Values.csv',dtype=float ,usecols=[2])
maximum = np.argmax(XRmax)
min_xmax = np.argmin(XRmax)
mins = min(XRmax)
#Find the largest and smallest R min vlaue
Rmins = np.loadtxt('Binary Log Values.csv', dtype=float, usecols=[3]) 
max_Rmin = np.argmax(Rmins)
maxR = max(Rmins)
min_Rmin = np.argmin(Rmins)
minRs = min(Rmins)



#For controlled sample (for the APOGEE random 713)
outfile2=open('Controlled Sample Log Values.lis', 'w')

Ax_max=[]
AR_min=[]
for k in range(len(ALX)):
    Ax_max.append(math.log10(ALX[k]))
    AR_min.append(math.log10(ALR[k]))
    outfile2.write(str(AmassID[k])+'\t'+str(Ax_max[k])+'\t'+str(AR_min[k])+ '\t' +str(ALX[k])+ '\t'+ str(ALR[k])+'\n')
outfile2.close()

perceBin = 0.05
binCount = int(perceBin*len(Ax_max))
gridz, Ax_range, Ay_range =np.histogram2d(Ax_max, AR_min, bins=binCount, range=None, normed=False, weights=None)

plt.figure(figsize=(8,6))
Z = plt.pcolor(Ax_range, Ay_range,gridz, cmap='RdBu')
plt.colorbar(label = 'Number of values')
plt.xlabel('log(max x-ranges)')
plt.ylabel('log(R min)')
plt.title('Controlled Sample of 713 stars from APOGEE-DR12')
plt.show()
plt.savefig('APOGEE_DR12 Candidate Log Values.jpg', bbox_inches='tight')
max_size  = np.argmax(gridz)
a1 = max_size%gridz.shape[0]
a2 = max_size/gridz.shape[0]

found_indi=find_index(Ax_max, AR_min, Ax_range, Ay_range, a2, a1)
print("Graham's algorithm found {0} indicies".format(found_indi))
print(AmassID[found_indi])
