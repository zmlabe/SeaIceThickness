"""
Script creates scatter plots for four quadrants of the Arctic between
PIOMAS and CryoSat-2

Author : Zachary Labe
Date : 8 August 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime
import scipy.stats as sts

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- PIOMAS Zoomed Regions (%s) ---' '\n' % titletime 

def piomasReader(directory,latmin,latmax,lonmin,lonmax,segment):
    """
    Reads piomas data for sea ice thickness over 1979-2015
    """
    
    ### Enter filename
    filename = 'piomas_regrid_March_19792015.nc'  
    
    ### Retrieve data
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
    
    ### Calculate lat/lon region
    xmask = (lat > latmin) & (lat < latmax)
    ymask = (lon > lonmin) & (lon < lonmax)
    
    mask = xmask[:] & ymask[:]
    latvals = np.where(mask == True)[0]
    lonvals = np.where(mask == True)[1]
    latvals = np.unique(latvals)
    lonvals = np.unique(lonvals)
    
    thk = sitp[:,latvals,:]
    thk = thk[:,:,lonvals]
    
    lat = lat[latvals,:]
    lat = lat[:,lonvals]
    lon = lon[latvals,:]
    lon = lon[:,lonvals]

    grid = '---> [[%s to %s N, %s to %s E]]' % (latmin,latmax,lonmin,lonmax)
    print 'Completed: PIOMAS data read (%s)!' % grid
    
    return lat,lon,thk
    
def satReader(directory,latmin,latmax,lonmin,lonmax):
    """
    Reads satellite data for joint CS-2/ICESat-J gridded (2004-2015)
    """
    
    ### Enter filename
    filename = 'satelliteJ_regrid_March_20042015.nc'
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    thkn = data.variables['thick'][:]
    data.close()
    
    ### Calculate lat/lon region
    xmask = (lat > latmin) & (lat < latmax)
    ymask = (lon > lonmin) & (lon < lonmax)
    
    mask = xmask[:] & ymask[:]
    latvals = np.where(mask == True)[0]
    lonvals = np.where(mask == True)[1]
    latvals = np.unique(latvals)
    lonvals = np.unique(lonvals)
    
    thk = thkn[:,latvals,:]
    thk = thk[:,:,lonvals]
    
    lat = lat[latvals,:]
    lat = lat[:,lonvals]
    lon = lon[latvals,:]
    lon = lon[:,lonvals]
    
    siti = thk[:6,:,:]
    sitc = thk[7:,:,:]

    grid = '---> [[%s to %s N, %s to %s E]]' % (latmin,latmax,lonmin,lonmax)
    print 'Completed: Satellite data read (%s)!' % grid
    
    return lat,lon,siti,sitc
    
### Call functions
latmin1 = 70
latmax1 = 89.9
lonmin1 = 20
lonmax1 = 90   
lat1,lon1,sitpi1 = piomasReader(directorydata,latmin1,latmax1,lonmin1,lonmax1,'icej')
lat1,lon1,sitpc1 = piomasReader(directorydata,latmin1,latmax1,lonmin1,lonmax1,'cryo')
lat1,lon1,siti1,sitc1 = satReader(directorydata,latmin1,latmax1,lonmin1,lonmax1)  

latmin2 = 70
latmax2 = 89.9
lonmin2 = 90
lonmax2 = 180 
lat2,lon2,sitpi2 = piomasReader(directorydata,latmin2,latmax2,lonmin2,lonmax2,'icej')
lat2,lon2,sitpc2 = piomasReader(directorydata,latmin2,latmax2,lonmin2,lonmax2,'cryo')
lat2,lon2,siti2,sitc2 = satReader(directorydata,latmin2,latmax2,lonmin2,lonmax2) 

latmin3 = 0
latmax3 = 89.9
lonmin3 = -180
lonmax3 = -90
lat3,lon3,sitpi3 = piomasReader(directorydata,latmin3,latmax3,lonmin3,lonmax3,'icej')
lat3,lon3,sitpc3 = piomasReader(directorydata,latmin3,latmax3,lonmin3,lonmax3,'cryo')
lat3,lon3,siti3,sitc3 = satReader(directorydata,latmin3,latmax3,lonmin3,lonmax3) 

latmin4 = 70
latmax4 = 89.9
lonmin4 = -90
lonmax4 = 0 
lat4,lon4,sitpi4 = piomasReader(directorydata,latmin4,latmax4,lonmin4,lonmax4,'icej')
lat4,lon4,sitpc4 = piomasReader(directorydata,latmin4,latmax4,lonmin4,lonmax4,'cryo')
lat4,lon4,siti4,sitc4 = satReader(directorydata,latmin4,latmax4,lonmin4,lonmax4) 

timex = np.arange(0,8,1)
timey = np.arange(0,8,1)

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

ax1 = plt.subplot(231)
plt.scatter(sitpi1[0,:,:],siti1[0,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[0,:,:],siti2[0,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[0,:,:],siti3[0,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[0,:,:],siti4[0,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
#ax1.plot(timex,line1,color='r',zorder=7)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'\textbf{2004}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax2 = plt.subplot(232)
plt.scatter(sitpi1[1,:,:],siti1[1,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[0,:,:],siti2[0,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[0,:,:],siti3[0,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[0,:,:],siti4[0,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'\textbf{2005}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax3 = plt.subplot(233)
plt.scatter(sitpi1[2,:,:],siti1[2,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[2,:,:],siti2[2,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[2,:,:],siti3[2,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[2,:,:],siti4[2,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{2006}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax4 = plt.subplot(234)
plt.scatter(sitpi1[3,:,:],siti1[3,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[3,:,:],siti2[3,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[3,:,:],siti3[3,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[3,:,:],siti4[3,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax4.set_xlim((0,7))
ax4.set_ylim((0,7))
ax4.set_xticklabels(map(str,np.arange(0,8,1)))
ax4.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax4, ['left', 'bottom'])
ax4.spines['top'].set_color('none')
ax4.spines['right'].set_color('none')
ax4.text(0,6.2,r'\textbf{2007}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax5 = plt.subplot(235)
plt.scatter(sitpi1[4,:,:],siti1[4,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[4,:,:],siti2[4,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[4,:,:],siti3[4,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[4,:,:],siti4[4,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax5.set_xlim((0,7))
ax5.set_ylim((0,7))
ax5.set_xticklabels(map(str,np.arange(0,8,1)))
ax5.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax5, ['left', 'bottom'])
ax5.spines['top'].set_color('none')
ax5.spines['right'].set_color('none')
ax5.text(0,6.2,r'\textbf{2008}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax6 = plt.subplot(236)
plt.scatter(sitpi1[5,:,:],siti1[5,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpi2[5,:,:],siti2[5,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpi3[5,:,:],siti3[5,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpi4[5,:,:],siti4[5,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax6.set_xlim((0,7))
ax6.set_ylim((0,7))
ax6.set_xticklabels(map(str,np.arange(0,8,1)))
ax6.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax6, ['left', 'bottom'])
ax6.spines['top'].set_color('none')
ax6.spines['right'].set_color('none')
ax6.text(0,6.2,r'\textbf{2009}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-9.4,-3.2,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
plt.text(-21,11.5,r'\textbf{sit( ICESat-J )}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'regionaltest.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################

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

ax1 = plt.subplot(231)
plt.scatter(sitpc1[0,:,:],sitc1[0,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpc2[0,:,:],sitc2[0,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpc3[0,:,:],sitc3[0,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpc4[0,:,:],sitc4[0,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'\textbf{2011}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax2 = plt.subplot(232)
plt.scatter(sitpc1[1,:,:],sitc1[1,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpc2[0,:,:],sitc2[0,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpc3[0,:,:],sitc3[0,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpc4[0,:,:],sitc4[0,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'\textbf{2012}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax3 = plt.subplot(233)
plt.scatter(sitpc1[2,:,:],sitc1[2,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpc2[2,:,:],sitc2[2,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpc3[2,:,:],sitc3[2,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpc4[2,:,:],sitc4[2,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{2013}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax4 = plt.subplot(234)
plt.scatter(sitpc1[3,:,:],sitc1[3,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpc2[3,:,:],sitc2[3,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpc3[3,:,:],sitc3[3,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpc4[3,:,:],sitc4[3,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax4.set_xlim((0,7))
ax4.set_ylim((0,7))
ax4.set_xticklabels(map(str,np.arange(0,8,1)))
ax4.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax4, ['left', 'bottom'])
ax4.spines['top'].set_color('none')
ax4.spines['right'].set_color('none')
ax4.text(0,6.2,r'\textbf{2014}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

ax5 = plt.subplot(235)
plt.scatter(sitpc1[4,:,:],sitc1[4,:,:],label='Region 1',color='seagreen',zorder=2,s=0.7)
plt.scatter(sitpc2[4,:,:],sitc2[4,:,:],label='Region 2',color='rosybrown',zorder=2,s=0.7)
plt.scatter(sitpc3[4,:,:],sitc3[4,:,:],label='Region 3',color='sandybrown',zorder=2,s=0.7)
plt.scatter(sitpc4[4,:,:],sitc4[4,:,:],label='Region 4',color='steelblue',zorder=2,s=0.7)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax5.set_xlim((0,7))
ax5.set_ylim((0,7))
ax5.set_xticklabels(map(str,np.arange(0,8,1)))
ax5.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax5, ['left', 'bottom'])
ax5.spines['top'].set_color('none')
ax5.spines['right'].set_color('none')
ax5.text(0,6.2,r'\textbf{2015}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.83,0.14),
                        frameon=False)

fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-.2,-3.2,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
plt.text(-12,11.5,r'\textbf{sit( CryoSat-2 )}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'regionaltest2.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Plot years
plt.figure()

region1 = sitpc1.copy()
region2 = sitpc2.copy()
region3 = sitpc3.copy()
region4 = sitpc4.copy()

region1[np.where(np.isnan(region1))] = 0.0
region1[np.where(np.isnan(region1))] = 0.0
region1[np.where(region1 != 0.0)] = 1.0
region1[np.where(region1 != 0.0)] = 1.0
region1[np.where(region1 != 1.0)] = np.nan
region1[np.where(region1 != 1.0)] = np.nan

region2[np.where(np.isnan(region2))] = 0.0
region2[np.where(np.isnan(region2))] = 0.0
region2[np.where(region2 != 0.0)] = 1.0
region2[np.where(region2 != 0.0)] = 1.0
region2[np.where(region2 != 1.0)] = np.nan
region2[np.where(region2 != 1.0)] = np.nan

region3[np.where(np.isnan(region3))] = 0.0
region3[np.where(np.isnan(region3))] = 0.0
region3[np.where(region3 != 0.0)] = 1.0
region3[np.where(region3 != 0.0)] = 1.0
region3[np.where(region3 != 1.0)] = np.nan
region3[np.where(region3 != 1.0)] = np.nan

region4[np.where(np.isnan(region4))] = 0.0
region4[np.where(np.isnan(region4))] = 0.0
region4[np.where(region4 != 0.0)] = 1.0
region4[np.where(region4 != 0.0)] = 1.0
region4[np.where(region4 != 1.0)] = np.nan
region4[np.where(region4 != 1.0)] = np.nan

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='darkgrey',linewidth=0.7)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.1,color='k',fontsize=9)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.1,color='k',fontsize=9)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream',alpha=1)

### Adjust maximum limits
values = np.arange(0,3,1)
info = 'Sea Ice Thickness'

### Plot filled contours    
cs = m.contourf(lon1,lat1,region1[0,:,:],
                values,latlon=True,colors='seagreen')
cs1 = m.contourf(lon2,lat2,region2[0,:,:],
                values,latlon=True,colors='rosybrown') 
cs2 = m.contourf(lon3,lat3,region3[0,:,:],
                values,latlon=True,colors='sandybrown')
cs3 = m.contourf(lon4,lat4,region4[0,:,:],
                values,latlon=True,colors='steelblue')  
                
plt.annotate(r'\textbf{I}',xy=(0.55,0.55),
             xycoords='axes fraction',fontsize=25,color='k')
plt.annotate(r'\textbf{II}',xy=(0.4,0.55),
             xycoords='axes fraction',fontsize=25,color='k')
plt.annotate(r'\textbf{III}',xy=(0.37,0.365),
             xycoords='axes fraction',fontsize=25,color='k')
plt.annotate(r'\textbf{IV}',xy=(0.55,0.365),
             xycoords='axes fraction',fontsize=25,color='k')

### Save figure
plt.savefig(directoryfigure +'quadrants.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Plot quadrants
### Calculate trends
varx1 = np.ravel(sitpi1)
vary1 = np.ravel(siti1)
mask1 = np.isfinite(varx1) & np.isfinite(vary1)
fit1 = np.polyfit(varx1[mask1],vary1[mask1],1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx1[mask1],
                                                          vary1[mask1])
m1 = fit1[0]
b1 = fit1[1]
linetest1 = m1*timex + b1

varx2 = np.ravel(sitpi2)
vary2 = np.ravel(siti2)
mask2 = np.isfinite(varx2) & np.isfinite(vary2)
fit2 = np.polyfit(varx2[mask2],vary2[mask2],1)
slope,intercept,r2,p_value,std_err = sts.stats.linregress(varx2[mask2],
                                                          vary2[mask2])
m2 = fit2[0]
b2 = fit2[1]
linetest2 = m2*timex + b2

varx3 = np.ravel(sitpi3)
vary3 = np.ravel(siti3)
mask3 = np.isfinite(varx3) & np.isfinite(vary3)
fit3 = np.polyfit(varx3[mask3],vary3[mask3],1)
slope,intercept,r3,p_value,std_err = sts.stats.linregress(varx3[mask3],
                                                          vary3[mask3])
m3 = fit1[0]
b3 = fit1[1]
linetest3 = m3*timex + b3

varx4 = np.ravel(sitpi4)
vary4 = np.ravel(siti4)
mask4 = np.isfinite(varx4) & np.isfinite(vary4)
fit4 = np.polyfit(varx4[mask4],vary4[mask4],1)
slope,intercept,r4,p_value,std_err = sts.stats.linregress(varx4[mask4],
                                                          vary4[mask4])
m4 = fit4[0]
b4 = fit4[1]
linetest4 = m4*timex + b4

fig = plt.figure()
ax1 = plt.subplot(221)
plt.scatter(sitpi1[0,:,:],siti1[0,:,:],label='2004',color='olive',zorder=2,s=1.5)
plt.scatter(sitpi1[1,:,:],siti1[1,:,:],label='2005',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpi1[2,:,:],siti1[2,:,:],label='2006',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpi1[3,:,:],siti1[3,:,:],label='2007',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpi1[4,:,:],siti1[4,:,:],label='2008',color='Aquamarine',zorder=2,s=1.5)
plt.scatter(sitpi1[5,:,:],siti1[5,:,:],label='2009',color='seagreen',zorder=2,s=1.5)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax1.plot(timex,linetest1,color='r',zorder=7)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'\textbf{Region I}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.23),
                        frameon=False)
ax1.text(2.95,0,r'r$^2$= %s' % round(r1,2),color='k',fontsize=7)

ax2 = plt.subplot(222)
plt.scatter(sitpi2[0,:,:],siti2[0,:,:],label='2004',color='olive',zorder=2,s=1.5)
plt.scatter(sitpi2[1,:,:],siti2[1,:,:],label='2005',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpi2[2,:,:],siti2[2,:,:],label='2006',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpi2[3,:,:],siti2[3,:,:],label='2007',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpi2[4,:,:],siti2[4,:,:],label='2008',color='Aquamarine',zorder=2,s=1.5)
plt.scatter(sitpi2[5,:,:],siti2[5,:,:],label='2009',color='seagreen',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.plot(timex,linetest2,color='r',zorder=7)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'\textbf{Region II}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.23),
                        frameon=False)
ax2.text(2.95,0,r'r$^2$= %s' % round(r2,2),color='k',fontsize=7)

ax3 = plt.subplot(223)
plt.scatter(sitpi3[0,:,:],siti3[0,:,:],label='2004',color='olive',zorder=2,s=1.5)
plt.scatter(sitpi3[1,:,:],siti3[1,:,:],label='2005',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpi3[2,:,:],siti3[2,:,:],label='2006',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpi3[3,:,:],siti3[3,:,:],label='2007',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpi3[4,:,:],siti3[4,:,:],label='2008',color='Aquamarine',zorder=2,s=1.5)
plt.scatter(sitpi3[5,:,:],siti3[5,:,:],label='2009',color='seagreen',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,linetest3,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{Region III}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.23),
                        frameon=False)
ax3.text(2.95,0,r'r$^2$= %s' % round(r3,2),color='k',fontsize=7)

ax3 = plt.subplot(224)
plt.scatter(sitpi4[0,:,:],siti4[0,:,:],label='2004',color='olive',zorder=2,s=1.5)
plt.scatter(sitpi4[1,:,:],siti4[1,:,:],label='2005',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpi4[2,:,:],siti4[2,:,:],label='2006',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpi4[3,:,:],siti4[3,:,:],label='2007',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpi4[4,:,:],siti4[4,:,:],label='2008',color='Aquamarine',zorder=2,s=1.5)
plt.scatter(sitpi4[5,:,:],siti4[5,:,:],label='2009',color='seagreen',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,linetest4,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{Region IV}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.23),
                        frameon=False)
ax3.text(2.95,0,r'r$^2$= %s' % round(r4,2),color='k',fontsize=7)

fig.subplots_adjust(wspace=.2)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-3,-3,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
plt.text(-10.3,10.5,r'\textbf{sit( ICESat-J )}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'regionalyrtest1.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Plot quadrants
varx1 = np.ravel(sitpc1)
vary1 = np.ravel(sitc1)
mask1 = np.isfinite(varx1) & np.isfinite(vary1)
fit1 = np.polyfit(varx1[mask1],vary1[mask1],1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx1[mask1],
                                                          vary1[mask1])
m1 = fit1[0]
b1 = fit1[1]
linetest1 = m1*timex + b1

varx2 = np.ravel(sitpc2)
vary2 = np.ravel(sitc2)
mask2 = np.isfinite(varx2) & np.isfinite(vary2)
fit2 = np.polyfit(varx2[mask2],vary2[mask2],1)
slope,intercept,r2,p_value,std_err = sts.stats.linregress(varx2[mask2],
                                                          vary2[mask2])
m2 = fit2[0]
b2 = fit2[1]
linetest2 = m2*timex + b2

varx3 = np.ravel(sitpc3)
vary3 = np.ravel(sitc3)
mask3 = np.isfinite(varx3) & np.isfinite(vary3)
fit3 = np.polyfit(varx3[mask3],vary3[mask3],1)
slope,intercept,r3,p_value,std_err = sts.stats.linregress(varx3[mask3],
                                                          vary3[mask3])
m3 = fit1[0]
b3 = fit1[1]
linetest3 = m3*timex + b3

varx4 = np.ravel(sitpc4)
vary4 = np.ravel(sitc4)
mask4 = np.isfinite(varx4) & np.isfinite(vary4)
fit4 = np.polyfit(varx4[mask4],vary4[mask4],1)
slope,intercept,r4,p_value,std_err = sts.stats.linregress(varx4[mask4],
                                                          vary4[mask4])
m4 = fit4[0]
b4 = fit4[1]
linetest4 = m4*timex + b4

fig = plt.figure()
ax1 = plt.subplot(221)
plt.scatter(sitpc1[0,:,:],sitc1[0,:,:],label='2011',color='olive',zorder=2,s=1.5)
plt.scatter(sitpc1[1,:,:],sitc1[1,:,:],label='2012',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpc1[2,:,:],sitc1[2,:,:],label='2013',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpc1[3,:,:],sitc1[3,:,:],label='2014',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpc1[4,:,:],sitc1[4,:,:],label='2015',color='Aquamarine',zorder=2,s=1.5)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax1.plot(timex,linetest1,color='r',zorder=7)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'\textbf{Region I}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.17),
                        frameon=False)
ax1.text(2.95,0,r'r$^2$= %s' % round(r1,2),color='k',fontsize=7)

ax2 = plt.subplot(222)
plt.scatter(sitpc2[0,:,:],sitc2[0,:,:],label='2011',color='olive',zorder=2,s=1.5)
plt.scatter(sitpc2[1,:,:],sitc2[1,:,:],label='2012',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpc2[2,:,:],sitc2[2,:,:],label='2013',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpc2[3,:,:],sitc2[3,:,:],label='2014',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpc2[4,:,:],sitc2[4,:,:],label='2015',color='Aquamarine',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.plot(timex,linetest2,color='r',zorder=7)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'\textbf{Region II}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.17),
                        frameon=False)
ax2.text(2.95,0,r'r$^2$= %s' % round(r2,2),color='k',fontsize=7)

ax3 = plt.subplot(223)
plt.scatter(sitpc3[0,:,:],sitc3[0,:,:],label='2011',color='olive',zorder=2,s=1.5)
plt.scatter(sitpc3[1,:,:],sitc3[1,:,:],label='2012',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpc3[2,:,:],sitc3[2,:,:],label='2013',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpc3[3,:,:],sitc3[3,:,:],label='2014',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpc3[4,:,:],sitc3[4,:,:],label='2015',color='Aquamarine',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,linetest3,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{Region III}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.17),
                        frameon=False)
ax3.text(2.95,0,r'r$^2$= %s' % round(r3,2),color='k',fontsize=7)

ax3 = plt.subplot(224)
plt.scatter(sitpc4[0,:,:],sitc4[0,:,:],label='2011',color='olive',zorder=2,s=1.5)
plt.scatter(sitpc4[1,:,:],sitc4[1,:,:],label='2012',color='rosybrown',zorder=2,s=1.5)
plt.scatter(sitpc4[2,:,:],sitc4[2,:,:],label='2013',color='sandybrown',zorder=2,s=1.5)
plt.scatter(sitpc4[3,:,:],sitc4[3,:,:],label='2014',color='steelblue',zorder=2,s=1.5)
plt.scatter(sitpc4[4,:,:],sitc4[4,:,:],label='2015',color='Aquamarine',zorder=2,s=1.5)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,linetest4,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'\textbf{Region IV}',color='k',fontsize=11)
plt.legend(shadow=False,fontsize=5,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.92,0.17),
                        frameon=False)
ax3.text(2.95,0,r'r$^2$= %s' % round(r4,2),color='k',fontsize=7)

fig.subplots_adjust(wspace=.2)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-3,-3,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
plt.text(-10.3,10.5,r'\textbf{sit( CryoSat-2 )}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'regionalyrtest2.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
### Plot regions by year for a single graph
### Region 1

def plotYrRegional(varx,vary,time,region):
    directoryfigure2 = '/home/zlabe/Desktop/RegionalYrScatter/%s/region%s/' % (time,region)
    
    if time == 'icesat':
        years = np.arange(2004,2010,1)
        labely = 'ICESat'
        color = 'seagreen'
        edgecolor = 'darkgreen'
    elif time == 'cryosat':
        years = np.arange(2011,2016,1)
        labely = 'CryoSat-2'
        color = 'steelblue'
        edgecolor = 'darkblue'
     
    timex = np.arange(0,8,1)
    timey = np.arange(0,8,1)     
    
    for i in xrange(len(varx)):
        vx= np.ravel(varx[i,:,:])
        vy = np.ravel(vary[i,:,:])
        mask = np.isfinite(vx) & np.isfinite(vy)
        fit = np.polyfit(vx[mask],vy[mask],1)
        slope,intercept,r,p_value,std_err = sts.stats.linregress(vx[mask],
                                                                 vy[mask])
        m = fit[0]
        b = fit[1]
        linetest = m*timex + b        
        
        fig = plt.figure()
        ax = plt.subplot(111)
        
        plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
        plt.scatter(varx[i,:,:],vary[i,:,:],label='%s' % years,
                    color=color,zorder=2,s=8,
                    edgecolor=edgecolor,linewidth=0.3)
        ax.plot(timex,linetest,color='r',zorder=3)
        
        ax.set_xlim((0,7))
        ax.set_ylim((0,7))
        ax.set_xticklabels(map(str,np.arange(0,8,1)))
        ax.set_yticklabels(map(str,np.arange(0,8,1)))
        adjust_spines(ax, ['left', 'bottom'])
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        
        plt.xlabel(r'sit( PIOMAS )',fontsize=13)
        plt.ylabel(r'sit( %s )' % labely,fontsize=13)
        plt.title(r'\textbf{Region %s -- %s}' % (region,years[i]),fontsize=13)
        ax.text(6.1,0.5,r'\textbf{r$^2$= %s}' % round(r,2),color='k',
                fontsize=11)
        
        fig.subplots_adjust(bottom=0.15)
                    
        plt.savefig(directoryfigure2 + '%s_regionScat_%s.png' % (time,years[i]),
                    dpi=300)
                    
plotYrRegional(sitpi1,siti1,'icesat',1)
plotYrRegional(sitpi2,siti2,'icesat',2)
plotYrRegional(sitpi3,siti3,'icesat',3)
plotYrRegional(sitpi4,siti4,'icesat',4)
plotYrRegional(sitpc1,sitc1,'cryosat',1)
plotYrRegional(sitpc2,sitc2,'cryosat',2)
plotYrRegional(sitpc3,sitc3,'cryosat',3)
plotYrRegional(sitpc4,sitc4,'cryosat',4)