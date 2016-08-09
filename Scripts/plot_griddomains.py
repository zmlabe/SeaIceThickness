"""
Script compares difference grid domains for CS-2, ICESat, and PIOMAS
 
Author : Zachary Labe
Date : 13 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime

### Define directories
directorydata1 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/Sat_thick/'  
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 2004
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
years = np.setdiff1d(years,[2010])      ### no satellite data in 2010
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- CS2/ICESat & PIOMAS Comparisons (%s) ---' '\n' % titletime 

def satReader(directory):
    """
    Reads satellite data for joint CS-2/ICESat gridded (2004-2015)
    """
    
    ### Enter filename
    filename = 'cs2icesat_regrid_mar_20042015.nc'    
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    thk_s = data.variables['thick'][:]
    data.close()
    
    print 'Completed: Satellite data read (March)!'
    
    return lats,lons,thk_s

def piomasReader(directory,month,years):
    """
    Reads piomas data for sea ice thickness over 1979-2015
    """
    
    ### Enter filename
    filename = 'piomas_regrid_sit_19792015.nc'   
    
    ### Month/Years extracted
    dateyr = now.year  
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    yearsp = np.arange(1979,2016)
    yearmin = years.min()
    yearmax = years.max()
    yearnone = 2010
    yearslice = np.where((yearsp <= yearmax) & (yearsp >= yearmin) & \
                (yearsp != yearnone))[0]
    
    ### Retrieve data
    data = Dataset(directory + filename)
    latp = data.variables['lat'][:]
    lonp = data.variables['lon'][:]
    thk_p = data.variables['newthickness'][yearslice,month,:,:]
    data.close()
    
    print 'Completed: PIOMAS data read (%s)!' % datemo
    
    return latp,lonp,thk_p
            
### Call functions
lats,lons,thk_s = satReader(directorydata2)  
latp,lonp,thk_p = piomasReader(directorydata1,2,years) 

thk_s[np.where(np.isnan(thk_s))] = 0.0
thk_p[np.where(np.isnan(thk_p))] = 0.0
thk_s[np.where(thk_s != 0.0)] = 1.0
thk_p[np.where(thk_p != 0.0)] = 1.0

thk_s[np.where(thk_s != 1.0)] = np.nan
thk_p[np.where(thk_p != 1.0)] = np.nan

### Define figure
fig = plt.figure()
ax = plt.subplot(111)

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.5)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.4,color='k',fontsize=9)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.4,color='k',fontsize=9)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

### Adjust maximum limits
values = np.arange(0,3,1)
info = 'Sea Ice Thickness'

### Plot filled contours    
cs = m.contourf(lons[:,:],lats[:,:],thk_p[0,:,:],
                values,latlon=True,colors='gold')
cs1 = m.contourf(lons,lats,thk_s[-5,:,:],
                values,latlon=True,colors='g',linewidths=1.3,alpha=0.5) 
cs2 = m.contourf(lons,lats,thk_s[2,:,:],
                values,latlon=True,colors='b',linewidths=1.3,alpha=0.5)                
    
fig.suptitle(r'\textbf{Grid Domains}, 100 km EASE',
         fontsize=15)
fig.subplots_adjust(top=0.85)

plt.annotate(r'\textbf{PIOMAS}',xy=(-0.5,1),
             xycoords='axes fraction',fontsize=22,color='gold')
             
plt.annotate(r'\textbf{CryoSat-2}',xy=(-0.5,0.85),
              xycoords='axes fraction',fontsize=22,color='green')
 
plt.annotate(r'\textbf{ICESat}',xy=(-0.5,0.7),
              xycoords='axes fraction',fontsize=22,color='blue')

### Save figure
plt.savefig(directoryfigure +'griddomains.png',dpi=300)

print 'Completed: Script done!'