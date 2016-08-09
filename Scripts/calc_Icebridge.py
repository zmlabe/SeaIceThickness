"""
Script reads ICESat-G data and regrids to EASE2.0 100 km grid. A netCDF
file is created.

Source : ftp://sidads.colorado.edu/pub/DATASETS/NSIDC0393_GLAS_SI_Freeboard_v01/
Author : Zachary Labe
Date : 20 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
from scipy.interpolate import griddata as g
import datetime
import numpy.ma as ma

### Define constants
years = np.arange(2009,2015)

### Define directories
directorydata = '/home/zlabe/Documents/SeaIceObs/Data/'
directorytest = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'IceBridge Data Read - %s' '\n' % titletime 

def readIcebridge(directory,years):
    """
    Function reads Icebridge data on EASE2.0 100 km for period of
    2009-2014. Spring data?
    """
    
    sit = np.empty((len(years),180,180))
    for i in xrange(len(years)):
        filename = 'icebridge_seaice_thickness_mean_spring_%s.nc' % years[i]
        
        data = Dataset(directory + filename)
        lats = data.variables['lat'][:]
        lons = data.variables['lon'][:]
        sit[i,:,:] = data.variables['thickness'][:]
        data.close()
        
    ### Check for missing data    
    sit[np.where(sit <= -9999.99)] = np.nan
    
    ### Read mean thickness
    filename2 = 'icebridge_seaice_thickness_mean_spring_2009to2014.nc'    
    
    data2 = Dataset(directory + filename2)
    meansit = data2.variables['thickness'][:]
    data2.close()
    
    ### Check for missing data
    meansit[np.where(meansit <= -9999.99)] = np.nan
    
    print 'Completed: Read Icebridge data!'
    return lats,lons,sit,meansit
    
def netcdfSIT(lats,lons,var,varmean):
    directory = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'
    name = 'icebridge_regrid_March_20092014.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Sea ice thickness processed by iceBridge onto' \
                         'EASE2.0 100km grids for spring 2009-2014.' \
                         'Mean thickness also included for period over' \
                         '2009-2014.'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('sit','f4',('years','lat','lon'))
    varnsmean = ncfile.createVariable('meansit','f4',('lat','lon'))
    
    ### Units
    varns.units = 'meters'
    varnsmean.units = 'meters'
    ncfile.title = 'IceBridge'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NSIDC, Andrew P. Barrett'
#    ncfile.references = ''
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    varnsmean[:] = varmean
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'

### Call functions
lats,lons,sit,meansit = readIcebridge(directorydata,years)
netcdfSIT(lats,lons,sit,meansit)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Practice plots to inspect regrids
def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
    
ice = meansit
lons = lons
lats = lats

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

### Define figure
fig = plt.figure()
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.25)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.25,color='k',fontsize=9)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.25,color='k',fontsize=9)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

### Adjust maximum limits
values = np.arange(0,7.1,0.1)  
ice[np.where(ice >= 7.0)] = 7.0
info = 'Sea Ice Thickness'

### Plot filled contours    
cs = m.contourf(lons[:,:],lats[:,:],ice[:,:],
                values,latlon=True,extend='max')
                  
### Set colormap                              
cs.set_cmap(colormapSIT())

### Set colorbar
cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)

cbar.set_ticks(np.arange(0,7.1,1))
cbar.set_ticklabels(map(str,np.arange(0,8,1)))    
cbar.set_label(r'\textbf{Thickness (meters)}')
    
fig.suptitle(r'\textbf{IceBridge (mean 2009-2014)}',
         fontsize=15)
fig.subplots_adjust(top=0.85)

### Save figure
plt.savefig(directorytest +'testicebridge.png',dpi=300)

print 'Completed: Script done!'