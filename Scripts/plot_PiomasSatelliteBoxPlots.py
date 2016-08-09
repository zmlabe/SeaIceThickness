"""
Script creates scatter plots between March data for SIT using PIOMAS
and ICESat-J and the submarine data
 
Author : Zachary Labe
Date : 27 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime
import scipy.stats as sts
import matplotlib.mlab as mlab
from matplotlib.ticker import NullFormatter

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
yearssub = np.arange(1986,1994+1,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- PIOMAS/Satellites Scatter Plots (%s) ---' '\n' % titletime 

def piomasReader(directory,segment,years):
    
    filename = 'piomas_regrid_March_19792015.nc'
    
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    
    if segment == 'sub':    # 1986-1994  
        timeslice = np.where((years >= 1986) & (years <= 1994))[0]
        sitp = data.variables['thick'][timeslice,:,:] 
    elif segment == 'icej':    # 2004-2009
        timeslice = np.where((years >= 2004) & (years <= 2009))[0]
        sitp = data.variables['thick'][timeslice,:,:]
    elif segment == 'cryo':    # 2011-2015
        timeslice = np.where((years >= 2011) & (years <= 2015))[0]
        sitp = data.variables['thick'][timeslice,:,:]
    else:
        sitp = data.variables['thick'][:,:,:]
        
    data.close()
    
    print 'Completed: PIOMAS data read!'
    return lat,lon,sitp
    
def subReader(directory):
    
    filename = 'sub_regrid_March_19861994.nc'    
    
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    sitb = data.variables['sit'][:]
    meansitb = data.variables['meansit'][:,:]
    data.close()
    
    print 'Completed: Submarine data read!'
    return lat,lon,sitb,meansitb
    
def icesatReader(directory):
    
    filename = 'satelliteJ_regrid_March_20042015.nc'
    
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    sit = data.variables['thick'][:]
    data.close()
    
    siti = sit[:6]
    sitc = sit[7:]
    
    print 'Completed: ICESat-J data read!'
    return lat,lon,siti,sitc
    
    
### Call functions
lat,lon,sitp = piomasReader(directorydata,'sub',years)
lat,lon,sitp2 = piomasReader(directorydata,'icej',years)
lat,lon,sitp3 = piomasReader(directorydata,'cryo',years)
lat,lon,sitb,meansitb = subReader(directorydata)
lat,lon,siti,sitc = icesatReader(directorydata)

varx1 = sitp2[0]
vary1 = siti[0]
mask1 = np.isfinite(varx1) & np.isfinite(vary1)
varx1 = varx1[mask1]
vary1 = vary1[mask1]

varx2 = sitp2[1]
vary2 = siti[1]
mask2 = np.isfinite(varx2) & np.isfinite(vary2)
varx2 = varx2[mask2]
vary2 = vary2[mask2]

varx3 = sitp2[2]
vary3 = siti[2]
mask3 = np.isfinite(varx3) & np.isfinite(vary3)
varx3 = varx3[mask3]
vary3 = vary3[mask3]

varx4 = sitp2[3]
vary4 = siti[3]
mask4 = np.isfinite(varx4) & np.isfinite(vary4)
varx4 = varx4[mask4]
vary4 = vary4[mask4]

varx5 = sitp2[4]
vary5 = siti[4]
mask5 = np.isfinite(varx5) & np.isfinite(vary5)
varx5 = varx5[mask5]
vary5 = vary5[mask5]

varx6 = sitp2[5]
vary6 = siti[5]
mask6 = np.isfinite(varx6) & np.isfinite(vary6)
varx6 = varx6[mask6]
vary6 = vary6[mask6]

dataq = [varx1,vary1,varx2,vary2,varx3,vary3,varx4,vary4,varx5,vary5,varx6,vary6]

timex = np.arange(1,13,1)
means = np.zeros((12))
maxes = np.zeros((12))
mins = np.zeros((12))
for i in xrange(len(means)):
    means[i] = np.nanmean(dataq[i])
    maxes[i] = np.percentile(dataq[i],75) + 1.5*(np.percentile(dataq[i],75) - np.percentile(dataq[i],25))
    mins[i] = np.percentile(dataq[i],25) - 1.5*(np.percentile(dataq[i],75) - np.percentile(dataq[i],25))
    
m,b = np.polyfit(timex,means,1)
mx,bax = np.polyfit(timex,maxes,1)
mn,bn = np.polyfit(timex,mins,1)

fig = plt.figure()
ax = plt.subplot(111)

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.tick_params(right='off',bottom='off')
ax.tick_params('both',length=5.5,width=2,which='major') 
ax.tick_params(axis='x',bottom='off',labelbottom='off') 

bx = plt.boxplot(dataq,0,'',patch_artist=True)
plt.plot(timex,m*timex + b,linestyle='-',color='k',linewidth=3,alpha=0.7)
plt.plot(timex,mx*timex + bax,linestyle='-',color='k',linewidth=3,alpha=0.7)
plt.plot(timex,mn*timex + bn,linestyle='-',color='k',linewidth=3,alpha=0.7)

plt.text(1.04,0,r'\textbf{2004}',fontsize=13)
plt.text(3.04,0,r'\textbf{2005}',fontsize=13)
plt.text(5.04,0,r'\textbf{2006}',fontsize=13)
plt.text(7.04,0,r'\textbf{2007}',fontsize=13)
plt.text(9.04,0,r'\textbf{2008}',fontsize=13)
plt.text(11.04,0,r'\textbf{2009}',fontsize=13)

plt.text(12.3,4.8,r'PIOMAS',fontsize=18,color='darkgrey',ha='right')
plt.text(12.3,4.4,r'ICESat-J',fontsize=18,color='steelblue',ha='right')

for i in bx['caps']:
    i.set(color='k',linewidth=0)
for whisker in bx['whiskers']:
    whisker.set(color='k',linestyle='-')
for box in bx['boxes'][::2]: # PIOMAS
    box.set(color='darkgrey')
for box in bx['boxes'][1::2]:
    box.set(color='steelblue')

plt.yticks(np.arange(0.5,5.5,0.5),map(str,np.arange(0.5,5.5,0.5)))    

plt.ylabel(r'Thickness (m)',fontsize=13)

plt.savefig(directoryfigure + 'boxtesta.png',dpi=300)