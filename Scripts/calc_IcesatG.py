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
years = np.arange(2003,2009)

### Define directories
directorygrid25 = '/home/zlabe/Surtsey/seaice_obs/sic/sic_rename/' 
directorygrid100 = '/home/zlabe/Documents/SeaIceObs/Data/'
directorydata = '/home/zlabe/Surtsey/seaice_obs/Icesat/icesatG/'
directorytest = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'ICESat-Goddard Read & Regrid - %s' '\n' % titletime 

def readIcesatG(directorygrid,directorydata,years):
    """
    Function reads ICESat data processed by NASA Goddard in a EASE2.0
    25 km grid
    """

    ### Read binary lat x lon arrays
    lons = np.fromfile(directorygrid + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directorygrid + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    sit = np.empty((len(years),136192))
    for i in xrange(len(years)):
        infile = directorydata + 'laser_thickness_mskd_%s.img' % years[i]
        with open(infile, 'rb') as fr:
            sit[i,:] = np.fromfile(fr, dtype='<f')
    
    ### Reshape array for [year x lat x lon]    
    sit = np.reshape(sit,(len(years),448,304))
    sit[np.where((sit == -4) | (sit == -3) | (sit == -2) | (sit == -1))] = np.nan
    
    print 'Completed: Read ICESat-G data!'
    return lats,lons,sit
    
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
    
    print 'Completed: Read grid data!' 
    return lats,lons
    
def regrid(lat1,lon1,lats,lons,var,years):
    """
    Interpolated on a 180x180 grid (latxlon) from CryoSat-2 using EASE2.0
    100 km grid
    """
    
    varn_re = np.reshape(var,(var.shape[0],(448*304)))       
    varn = np.empty((var.shape[0],lats.shape[0],lons.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(varn.shape[0]):
        z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,:],
              (lats,lons),method='linear')
        varn[i,:,:] = z
        print 'Completed: Year %s Regridding---' % (years[i])
    print 'Completed: Done regridding process!'
    return varn
    
def netcdfSIT(lats,lons,var):
    directory = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'
    name = 'icesatG_regrid_March_20032008.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Sea ice thickness processed by NASA-G and now' \
                         'regridded on an EASE2.0 100 km grid for the' \
                         'period of March 2003-2008'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('sit','f4',('years','lat','lon'))
    
    ### Units
    varns.units = 'meters'
    ncfile.title = 'ICESat-G'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NASA-G'
    ncfile.references = 'Donghui Yi, H. Zwally'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'

### Call functions
#lats,lons,sit = readIcesatG(directorygrid25,directorydata,years)
#latsnew,lonsnew = gridData(directorygrid100)
#sitn = regrid(lats,lons,latsnew,lonsnew,sit,years)
#netcdfSIT(latsnew,lonsnew,sitn)

print 'Completed: Script done!'

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Practice plots to inspect regrids

#def colormapSIT():
#    cmap1 = plt.get_cmap('BuPu')
#    cmap2 = plt.get_cmap('RdPu_r')
#    cmap3 = plt.get_cmap('gist_heat_r')
#    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
#    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
#    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
#    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
#    return cms_sit
#    
#ice = sit[0,:,:]
#lons = lons
#lats = lats
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
#values = np.arange(0,7.1,0.1)  
#ice[np.where(ice >= 7.0)] = 7.0
#info = 'Sea Ice Thickness'
#
#### Plot filled contours    
#cs = m.contourf(lons[:,:],lats[:,:],ice[:,:],
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
#fig.suptitle(r'\textbf{[ICESat-G]}',
#         fontsize=15)
#fig.subplots_adjust(top=0.85)
#
#### Save figure
#plt.savefig(directorytest +'testice.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Test difference between ICESat-G and ICESat-J

def readICESatJ():
    
    directory = '/home/zlabe/Documents/SeaIceObs/Data/'
    years = np.arange(2004,2010)
    
    thickness_icesat = np.empty((len(years),180,180))
    for i in xrange(len(years)):
        files = 'icesat_seaice_thickness_mean_spring_%s.nc' % years[i]
        filename = directory + files
        data = Dataset(filename)
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        thickness_icesat[i,:,:] = data.variables['thick'][:]
        data.close()
        
    thickness_icesat[np.where(thickness_icesat == 9999.)] = np.nan
    lon[np.where(lon == -999.)] = np.nan
    lat[np.where(lat == -999.)] = np.nan
    
    print 'Completed: ICESat-J Data Read!'
        
    return thickness_icesat,lat,lon
    
sitj,lat,lon = readICESatJ()   

add = np.empty((1,180,180))
add.fill(np.nan)

sit_j = np.append(add,sitj,axis=0)
sit_g = np.append(sitn,add,axis=0)

def readPiomas(directory):
    """
    Reads netCDF4 PIOMAS data
    """

    files = 'piomas_regrid_sit_19792015.nc'        
    
    data = Dataset(directory + files)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    sit = data.variables['newthickness'][:,2,:,:]
    data.close()

    print 'Completed: PIOMAS Read data!'  
    
    return lats,lons,sit
lats,lons,sitp = readPiomas('/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/')

thickness = sit_j - sit_g
thickness2 = sit_g - sitp[-13:-6]
thickness3 = sit_j - sitp[-13:-6]
lons = lons
lats = lats

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

yearnew = np.arange(2003,2010,1)

### Define figure
fig = plt.figure()    
for i in xrange(len(thickness)):
    ax = plt.subplot(2,4,i+1)    
    m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
    m.drawmapboundary(fill_color = 'white')
    m.drawcoastlines(color = 'darkgrey',linewidth=0.2)
    m.drawlsmask(land_color='darkgrey',ocean_color='snow')
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.25)
    m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.25)
    
    cs = m.contourf(lon,lat,thickness3[i,:,:],np.arange(-3,3.1,0.25),
                    latlon=True,extend='both')
    cs.set_cmap('seismic_r')
    ax.text(0.89,0.95,r'\textbf{%s}' % (yearnew[i]),size='9',
            horizontalalignment='center',backgroundcolor='w',
            verticalalignment='center',bbox=dict(facecolor='w',
            edgecolor='k',alpha=0.9),transform=ax.transAxes)
            
plt.tight_layout()
fig.subplots_adjust(bottom=0.13)
cbar_ax = fig.add_axes([0.15,0.1,0.7,0.04])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac='auto',drawedges=True)

cbar.set_label(r'Thickness (meters)')
cbar.set_ticks(np.arange(-3,4,1))
cbar.set_ticklabels(map(str,np.arange(-3,4,1)))                  
fig.suptitle(r'\textbf{[ICESat-J -- PIOMAS]}',
             fontsize=20)
fig.subplots_adjust(top=0.93)
fig.subplots_adjust(wspace=0.02)
fig.subplots_adjust(hspace=-0.09)
    
plt.savefig(directorytest +'diffe3.png',dpi=300)