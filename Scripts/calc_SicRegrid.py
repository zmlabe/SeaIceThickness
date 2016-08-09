"""
Script reads Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS
Passive Microwave Data, Version 1 binary files for select variables and 
regrids according to selected grid style (e.g., NSIDC EASE grid data). 
 
Source : https://nsidc.org/data/nsidc-0051#
Author : Zachary Labe
Date : 14 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
from scipy.interpolate import griddata as g
import datetime
#import cartopy.crs as ccrs
import numpy.ma as ma

### Define directories
directorygrid = '/home/zlabe/Documents/SeaIceObs/Data/'
directorydata = '/home/zlabe/Surtsey/seaice_obs/sic/sic_rename/'    
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1993               # first time includes 12 months
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = np.arange(1,13,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'Satellite SIC Read & Regrid - %s' '\n' % titletime 

def SicRead(directory,years):
    """
    Reads binary sic data
    """
    
    ### Read binary lat x lon arrays
    lons = np.fromfile(directory + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directory + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    ### Read binary sea ice concentration
    ice = np.empty((len(years),len(months),136192))
    for i in xrange(len(years)):
        for j in xrange(len(months)):
            
            if months[j] < 10:
                filename = 'sic_%s_0%s.bin' % (years[i],months[j])
            else:
                filename = 'sic_%s_%s.bin' % (years[i],months[j])
            infile = directory + filename
            
            with open(infile, 'rb') as fr:
                hdr = fr.read(300)
                ice[i,j,:] = np.fromfile(fr, dtype=np.uint8)
    
    ice = np.reshape(ice,(ice.shape[0],ice.shape[1],448,304))
    ice = ice/250               # Scale Factor
    
    ### Assign mask
    mask = np.fromfile(directory + 'gsfc_25n.msk',dtype='int8')
    mask = np.reshape(mask,(448,304))
    mask = mask.astype(float)
    mask[np.where(mask == 1.0)] = np.nan
    mask[np.where(mask == 0.0)] = 1.0
    ice = ice*mask
    
    print 'Completed: Read SIC data!' 
    return lats,lons,ice,hdr,mask

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
    
def regrid(lat1,lon1,lats,lons,var,years):
    """
    Interpolated on a 180x180 grid (latxlon) from CryoSat-2 using EASE2.0
    100 km grid
    """
    
    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(448*304)))       
    varn = np.empty((var.shape[0],var.shape[1],lats.shape[0],lons.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(varn.shape[0]):
        for j in xrange(varn.shape[1]):
            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],
                  (lats,lons),method='linear')
            varn[i,j,:,:] = z
        print 'Completed: Year %s Regridding---' % (years[i])
    print 'Completed: Done regridding process!'
    return varn

def netcdfSIC(lats,lons,var):
    directory = '/home/zlabe/Surtsey/seaice_obs/sic/'
    name = 'nsidc_regrid_sic_19932015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Sea Ice Concentrations from Nimbus-7 SMMR and' \
                         'DMSP SSM/I-SSMIS Passive Microwave Data,' \
                         'Version 1 -- files regridded with EASE25'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('months',var.shape[1])
    ncfile.createDimension('lat',var.shape[2])
    ncfile.createDimension('lon',var.shape[3])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('sic','f4',('years','months','lat','lon'))
    
    ### Units
    varns.units = 'fraction (%)'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    months[:] = list(xrange(var.shape[1]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'
    
### Call functions
lats1,lons1,ice,hdr,mask = SicRead(directorydata,years)
#lats,lons = gridData(directorygrid)
#sic = regrid(lats1,lons1,lats,lons,ice,years)
#netcdfSIC(lats,lons,sic)

print 'Completed: Script finished!'

###########################################################################
###########################################################################
###########################################################################
#def colormapSIC():
#    cmap = plt.get_cmap('RdPu')
#    cmaplist = [cmap(i) for i in xrange(0,cmap.N-20)]
#    cms_sic = c.ListedColormap(cmaplist)
#    return cms_sic
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
#m.drawcoastlines(color='k',linewidth=0.5)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.5,color='k',fontsize=9)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0.5,color='k',fontsize=9)
#m.drawlsmask(land_color='peru',ocean_color='mintcream')
#
#### Adjust maximum limits
#values = np.arange(0.2,1.05,0.05)
#sic[np.where(sic >= 1.0)] = 1.0
#info = 'SIC'
#
#### Plot filled contours    
#cs = m.contourf(lons[:,:],lats[:,:],sic[0,0,:,:],
#                values,latlon=True,extend='max')
#                  
#### Set colormap                              
#cs.set_cmap(colormapSIC())
#
#### Set colorbar
#cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)
#
#cbar.set_ticks(np.arange(0.2,1.1,.1))
#cbar.set_ticklabels(map(str,np.arange(0.2,1.1,0.1)))
#cbar.set_label(r'\textbf{SIC Fraction}')
#    
#fig.suptitle(r'\textbf{SIC Test}',
#         fontsize=15)
#fig.subplots_adjust(top=0.85)
#
#### Save figure
#plt.savefig(directoryfigure +'testsic.png',dpi=300)