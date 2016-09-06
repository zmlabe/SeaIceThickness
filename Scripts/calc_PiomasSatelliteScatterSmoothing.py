"""
Script creates scatter plots between March data for SIT using PIOMAS
and ICESat-J and CryoSat-2. Subplots are for each individual year 
over the 2004-2015 time frame. Script also looks at techniques for smoothing
to understand the linear fit.
 
Author : Zachary Labe
Date : 15 August 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import statsmodels.api as sm

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
directoryfigure = '/home/zlabe/Desktop/Smoothing/'
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

print '\n' '--- PIOMAS/Satellites Yearly Scatter Plots & Smoothing (%s) ---' '\n' % titletime 

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
lat,lon,sitpi = piomasReader(directorydata,'icej',years)
lat,lon,sitpc = piomasReader(directorydata,'cryo',years)
lat,lon,siti,sitc = icesatReader(directorydata)

timex = np.arange(0,8,1)
timey = np.arange(0,8,1)

### Smoothing test
varx = np.ravel(sitpi)
vary = np.ravel(siti)
mask = np.isfinite(varx) & np.isfinite(vary)

varx = varx[mask]
vary = vary[mask]

test = np.arange(0,7,7./len(varx))

smoothed = sm.nonparametric.lowess(vary,varx,it=0,frac=0.5)

plt.scatter(varx,vary)
plt.plot(smoothed[:,0],smoothed[:,1],color='r')

### ICESat-J yearly plots
fig = plt.figure()

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

def plotYrScatter(varix,variy,time):
    if time == 'icesat':
        years = np.arange(2004,2010,1)
        labely = 'ICESat-J'
        color = 'seagreen'
        edgecolor = 'darkgreen'
        cordx1 = -9.4
        cordy1 = -3.2
        cordx2 = -21
        cordy2 = 11.5
    elif time == 'cryosat':
        years = np.arange(2011,2016,1)
        labely = 'CryoSat-2'
        color = 'steelblue'
        edgecolor = 'darkblue'
        cordx1 = -0.2
        cordy1 = -3.2
        cordx2 = -12
        cordy2 = 11.5
    
    fig = plt.figure()
    for i in xrange(len(varix)):
        ax = fig.add_subplot(2,3,i+1)
        
        ### Calculate smoothing
        varx = np.ravel(varix[i,:,:])
        vary = np.ravel(variy[i,:,:])
        mask = np.isfinite(varx) & np.isfinite(vary)  
        varx = varx[mask]
        vary = vary[mask]
        
        smoothed = sm.nonparametric.lowess(vary,varx,it=0,frac=0.5)
        
        plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
        plt.scatter(varx,vary,s=10,color=color,edgecolor=edgecolor,
                    linewidth=0.3,zorder=2)
        plt.plot(smoothed[:,0],smoothed[:,1],color='r',zorder=3)
        
        ax.set_xlim((0,7))
        ax.set_ylim((0,7))
        ax.set_xticklabels(map(str,np.arange(0,8,1)))
        ax.set_yticklabels(map(str,np.arange(0,8,1)))
        adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        
        ax.text(0,6.2,r'%s' % years[i],color=color,fontsize=11)
     
    fig.subplots_adjust(wspace=.3)
    fig.subplots_adjust(hspace=.4)
    fig.subplots_adjust(bottom=0.15) 
    
    
    plt.text(cordx1,cordy1,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
    plt.text(cordx2,cordy2,r'\textbf{sit( %s )}' % labely,fontsize=16,rotation=90)
    plt.savefig(directoryfigure + 'yrsmooth_%s.png' % time,dpi=300)
    
plotYrScatter(sitpi,siti,'icesat')
plotYrScatter(sitpc,sitc,'cryosat')