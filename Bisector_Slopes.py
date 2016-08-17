import numpy as np
from matplotlib import pyplot as plt
import math
import matplotlib.mlab as mlab
from astropy.io import fits
import pylab as pl
import matplotlib.lines as mlines

# Reading in the txt associated with the R and x-range values 
SB_x=np.loadtxt('SB9_Bisector_pts.csv', skiprows=1, dtype =str, usecols=[0,1,2,3])
DrewsList_x=np.loadtxt('DrewsList_Bisector_pts.lis', skiprows=1, dtype=str, usecols=[0,1,2,3])
DrewsList2_x=np.loadtxt('DrewsList2_Bisector_pts.lis', skiprows=1, dtype=str, usecols=[0,1,2,3])
DR12_random_x=np.loadtxt('APOGEE_Random_Bisector_pts.lis',skiprows=1, dtype=str, usecols=[0,1,2,3] )

twoMASSID=np.loadtxt('DrewsList_Bisector_pts.lis',skiprows=1,dtype = str, usecols=[0])
MassID = np.loadtxt('SB9_Bisector_pts.csv', skiprows=1, dtype=str, usecols=[0])
print(len(twoMASSID))
y = [0.1,0.3,0.5]

SB1 = SB_x[0:3]
SB2 = SB_x[3]
SB3 = SB_x[4:13]
SB4 = SB_x[13:21]
SB5 = SB_x[21:30]
SB6 = SB_x[30:39]
SB7 = SB_x[39:43]
SB8 = SB_x[43:47]
SB9 = SB_x[47:50]
SB10 = SB_x[50:53]
SB11 = SB_x[53:56]
SB12 = SB_x[56:58]
SB12 = SB_x[58:61]
SB13 = SB_x[61:66]
SB14 = SB_x[66:72]
SB15 = SB_x[72]
SB16 = SB_x[73:77]
SB17 = SB_x[77:88]
SB18 = SB_x[88:91]
SB19 = SB_x[91:94]
SB20 = SB_x[94:97]
SB21 = SB_x[97:100]
SB22 = SB_x[100:103]
SB23 = SB_x[103:106]
SB24 = SB_x[106:109]
SB25 = SB_x[109:112]
SB26 = SB_x[112:115]
SB27 = SB_x[115:127]
SB28 = SB_x[127:149]
SB29 = SB_x[149:152]
SB30 = SB_x[152:160]
SB31 = SB_x[160:168]
SB32 = SB_x[168:171]
SB33 = SB_x[171:179]
SB34 = SB_x[179:182]
SB35 = SB_x[182:185]
SB36 = SB_x[185:188]
SB37 = SB_x[188:191]
SB38 = SB_x[191:195]
SB39 = SB_x[212:216]
SB40 = SB_x[217:220]
SB41 = SB_x[221:224]
SB42 = SB_x[221:234]
SB43 = SB_x[234:238]
SB44 = SB_x[238:241]
SB45 = SB_x[241:244]
SB46 = SB_x[244:247]
SB47 = SB_x[247:251]
SB48 = SB_x[251:254]
SB49 = SB_x[254:257]
SB50 = SB_x[257:262]
SB51 = SB_x[262:275]
SB52 = SB_x[275:288]
SB53 = SB_x[288:291]
SB54 = SB_x[291:295]
SB55 = SB_x[295:313]
SB56 = SB_x[313:316]
SB57 = SB_x[316:320]
SB58 = SB_x[320:325]
SB59 = SB_x[325:330]
SB60 = SB_x[330:332]
SB61 = SB_x[332:340]
SB62 = SB_x[340:343]
SB63 = SB_x[343:356]
SB64 = SB_x[356:358]
SB65 = SB_x[358]
SB66 = SB_x[359:362]
SB67 = SB_x[362:365]
SB68 = SB_x[365:368]
SB69 = SB_x[368:271]
SB70 = SB_x[371:374]
SB71 = SB_x[374:377]
SB72 = SB_x[377:380]

#Plotting bisectors of only three points
plt.figure(figsize=(8,8))
slope=[]
poly_coeff=[]
for i in range(len(SB_x)):
    
    plt.figure(figsize=(8,8))
    plt.plot(SB_x[i],y,'ro')
    plt.xlabel('CCF lag')
    plt.ylabel('Agreement')
    plt.title('Bisectors for '+str(MassID[i]))
    plt.plot(SB_x[i], np.poly1d(np.polyfit(SB_x[i], y, 1))(SB_x[i]));
    slope.append(np.polyfit(SB_x[i],y,1));
    m = np.polyfit(SB_x[i],y,1);

    plt.plot(SB_x[i],np.poly1d(np.polyfit(SB_x[i],y,4))(SB_x[i]));
    poly_coeff.append(np.poly1d(np.polyfit(SB_x[i],y,4)));
    
    
#Writing a csv to hold the 2MASSID and bisector linear coefficients for three points
outfile2 = open('SB9_Coefficients_for_3pts.csv','w')
for i in range(len(SB_x)):
    outfile2.write(str(MassID[i])+'\t'+str(slope[i][0])+'\n')

outfile2.close()

