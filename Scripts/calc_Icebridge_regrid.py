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
directorygrid = '/home/zlabe/Documents/SeaIceObs/Data/'
directorydata = '/home/zlabe/Surtsey/seaice_obs/Icebridge/mission2015/'
directorytest = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'IceBridge Data Read and Regrid - %s' '\n' % titletime 

def gridData(directory):
    """
    Reads CryoSat-2 data that is AWI processed onto EASE100 grid
    """
    
    files = 'cryoSat_seaice_thickness_month_2011.nc' 
    filename = directory + files
    data1 = Dataset(filename)
    lats = data1.variables['lat'][:]
    lons = data1.variables['lon'][:]
    data1.close()
    lons[np.where(lons == -999)] = np.nan
    lats[np.where(lats == -999)] = np.nan
    
    print 'Completed: Read grid data!' '\n'
    return lats,lons

lats,lons = gridData(directorygrid)

filename = 'OIB_20150326_IDCSI2.txt'
filename2 = 'OIB_20150327_IDCSI2.txt'
filename3 = 'OIB_20150401_IDCSI2.txt'
filename4 = 'OIB_20150329_IDCSI2.txt'
lat1,lon1,sit = np.genfromtxt(directorydata + filename,delimiter=',',
                            unpack=True,skip_header=1,usecols=[0,1,2])
lat2,lon2,sit2 = np.genfromtxt(directorydata + filename2,delimiter=',',
                            unpack=True,skip_header=1,usecols=[0,1,2])  
lat3,lon3,sit3 = np.genfromtxt(directorydata + filename3,delimiter=',',
                            unpack=True,skip_header=1,usecols=[0,1,2])
lat4,lon4,sit4 = np.genfromtxt(directorydata + filename4,delimiter=',',
                            unpack=True,skip_header=1,usecols=[0,1,2])                             
sit[np.where(sit <= -9999)] = np.nan   
sit2[np.where(sit2 <= -9999)] = np.nan
sit3[np.where(sit3 <= -9999)] = np.nan
sit4[np.where(sit4 <= -9999)] = np.nan

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)       

x,y = m(lon1,lat1)  
x2,y2 = m(lon2,lat2)     
x3,y3 = m(lon3,lat3)   
x4,y4 = m(lon4,lat4)          
                            
#z = g((lat1,lon1),sit,(lats,lons),method='nearest')

    
#def netcdfSIT(lats,lons,var,varmean):
#    directory = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'
#    name = 'icebridge_regrid_March_20092014.nc'
#    filename = directory + name
#    ncfile = Dataset(filename,'w',format='NETCDF4')
#    ncfile.description = 'Sea ice thickness processed by iceBridge onto' \
#                         'EASE2.0 100km grids for spring 2009-2014.' \
#                         'Mean thickness also included for period over' \
#                         '2009-2014.'
#    
#    ### Dimensions
#    ncfile.createDimension('years',var.shape[0])
#    ncfile.createDimension('lat',var.shape[1])
#    ncfile.createDimension('lon',var.shape[2])
#    
#    ### Variables
#    years = ncfile.createVariable('years','f4',('years'))
#    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
#    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
#    varns = ncfile.createVariable('sit','f4',('years','lat','lon'))
#    varnsmean = ncfile.createVariable('meansit','f4',('lat','lon'))
#    
#    ### Units
#    varns.units = 'meters'
#    varnsmean.units = 'meters'
#    ncfile.title = 'IceBridge'
#    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
#    ncfile.source = 'NSIDC, Andrew P. Barrett'
##    ncfile.references = ''
#    
#    ### Data
#    years[:] = list(xrange(var.shape[0]))
#    latitude[:] = lats
#    longitude[:] = lons
#    varns[:] = var
#    varnsmean[:] = varmean
#    
#    ncfile.close()
#    print 'Completed: Created netCDF4 File!'
#
#### Call functions
##netcdfSIT(lats,lons,sit,meansit)
#
############################################################################
############################################################################
############################################################################
############################################################################
### Practice plots to inspect regrids
#ice = z
#
#
def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
    
fig = plt.figure()
ax = plt.subplot(111)
plt.hexbin(x2,y2,C=sit2,gridsize=100,cmap=colormapSIT(),rasterized=True)
plt.hexbin(x3,y3,C=sit3,gridsize=100,cmap=colormapSIT(),rasterized=True)
plt.hexbin(x4,y4,C=sit4,gridsize=100,cmap=colormapSIT(),rasterized=True)
plt.hexbin(x,y,C=sit,gridsize=100,cmap=colormapSIT(),rasterized=True)
plt.colorbar()
plt.title(r'IceBridge March 2015')

#
#### Call parameters
#plt.rcParams['text.usetex']=True
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Avant Garde'
#
#### Define figure
#fig = plt.figure()
#ax = plt.subplot(111)
#
#m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
#    
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='k',linewidth=0.25)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.25,color='k',fontsize=9)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0.25,color='k',fontsize=9)
#m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#
#### Adjust maximum limits
#values = np.arange(0,7.1,0.1)  
#ice[np.where(ice >= 7.0)] = 7.0
#info = 'Sea Ice Thickness'
#
#### Plot filled contours    
#cs = m.contourf(lons,lats,ice,
#                values,latlon=True,extend='max')
#                  
#### Set colormap                              
#cs.set_cmap(colormapSIT())
#
#### Set colorbar
#cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)
#
#cbar.set_ticks(np.arange(0,7.1,1))
#cbar.set_ticklabels(map(str,np.arange(0,8,1)))    
#cbar.set_label(r'\textbf{Thickness (meters)}')
#    
#fig.suptitle(r'\textbf{IceBridge (mean 2009-2014)}',
#         fontsize=15)
#fig.subplots_adjust(top=0.85)
#
#### Save figure
plt.savefig(directorytest +'testicebridge2.png',dpi=300)

print 'Completed: Script done!'