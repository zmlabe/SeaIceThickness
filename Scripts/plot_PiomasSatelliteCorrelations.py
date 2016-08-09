"""
Script calculates pearson temporal correlations between PIOMAS and
submarine/ICESat-J/CryoSat-2 data
 
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
import iris as ir
import iris.quickplot as qplt

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

print '\n' '--- PIOMAS/Satellites Correlation Plots (%s) ---' '\n' % titletime 

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
lat,lon,sitpb = piomasReader(directorydata,'sub',years)
lat,lon,sitpi = piomasReader(directorydata,'icej',years)
lat,lon,sitpc = piomasReader(directorydata,'cryo',years)
lat,lon,sitb,meansitb = subReader(directorydata)
lat,lon,siti,sitc = icesatReader(directorydata)

def Corr(x,y):
    """
    Calculates pearson correlation coefficent between two variables
    
    
    Parameters
    ----------
    x : dependent variable
    y : independent variable
    
    Returns
    ----------
    cocoeff1 : pearson correlation coefficient 
    cocoeff2 : not important correlation coefficient 
    """
    
    cocoeff1 = np.empty((y.shape[1],y.shape[2]))
    cocoeff2 = np.empty((y.shape[1],y.shape[2]))
    for i in xrange(y.shape[1]):
        for j in xrange(y.shape[2]):
            cocoeff1[i,j],cocoeff2[i,j] = sts.pearsonr(x[:,i,j],y[:,i,j])
    
    print 'Completed: Correlation calculations!'
            
    return cocoeff1, cocoeff2 
    
corr1_b,corr2_b = Corr(sitpb,sitb)
corr1_i,corr2_i = Corr(sitpi,siti)
corr1_c,corr2_c = Corr(sitpc,sitc)

### Plot Correlations
fig = plt.figure()

cmap = plt.cm.get_cmap('brewer_RdBu_11')  

plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde' 

ax = fig.add_subplot(1,3,1)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
cs1 = m.contourf(lon,lat,corr1_b,np.arange(-1,1.1,0.1),latlon=True)
cs1.set_cmap(cmap)

ax = fig.add_subplot(1,3,2)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
cs1 = m.contourf(lon,lat,corr1_i,np.arange(-1,1.1,0.1),latlon=True)
cs1.set_cmap(cmap)

ax = fig.add_subplot(1,3,3)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.3,color='k',fontsize=3)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
cs1 = m.contourf(lon,lat,corr1_c,np.arange(-1,1.1,0.1),latlon=True)
cs1.set_cmap(cmap)

cbar_ax = fig.add_axes([0.26,0.2,0.5,0.045])                
cbar = fig.colorbar(cs1,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)
cbar.set_ticks(np.arange(-1,1.2,0.2))
cbar.set_ticklabels([r'-1.0',r'-0.8',r'-0.6',r'-0.4',r'-0.2',r'0.0',
                   r'0.2',r'0.4',r'0.6',r'0.8',r'1.0'])
cbar.set_label(r'\textbf{Correlation Coefficient}')

plt.annotate(r'\textbf{Submarine}',xy=(-0.23,11.5),xycoords='axes fraction',
     fontsize=13)
plt.annotate(r'\textbf{ICESat-J}',xy=(0.377,11.5),xycoords='axes fraction',
     fontsize=13)
plt.annotate(r'\textbf{CryoSat-2}',xy=(0.91,11.5),xycoords='axes fraction',
     fontsize=13)

fig.subplots_adjust(wspace=0.3)

plt.savefig(directoryfigure + 'corrtest.png',dpi=300)