"""
Script creates time series for month of March using satellite data and
modeled data from PIOMAS for sea ice thickness
 
Source 1 : ftp://sidads.colorado.edu/pub/projects/SIPN/seaice_thickness/
Source 2 : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
Author : Zachary Labe
Date : 19 July 2016
"""
### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata as g
import datetime
import numpy.ma as ma

### Define directories
directorygrid = '/home/zlabe/Documents/SeaIceObs/Data/'
directorydata_p = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'    
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 
directorynew = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'March SIT Data Set - %s' '\n' % titletime 

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
    
def netcdfPiomas(lats,lons,var,directory):
    name = 'piomas_regrid_March_19792015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'PIOMAS Sea ice thickness reanalysis from 1979-2015 ' \
                        'interpolated on a 180x180 grid (latxlon)' \
                        'of NSIDC EASE100' 
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('thick','f4',('years','lat','lon'))
    
    ### Metrics
    varns.units = 'meters'
    ncfile.title = 'PIOMAS March SIT'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'University of Washington'
    ncfile.references = '[Zhang and Rothrock, 2003]'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'
    
### Call functions
lats_p,lons_p,sit_p = readPiomas(directorydata_p)
#netcdfPiomas(lats_p,lons_p,sit_p,directorynew)

###########################################################################
###########################################################################
###########################################################################
###########################################################################

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
    
def readCryoSat():

    directory = '/home/zlabe/Documents/SeaIceObs/Data/'
    
    filename = directory + 'cryosat_March_20112015.nc'
    data = Dataset(filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    thickness_cryosat = data.variables['newthickness'][:]
    data.close()
    
    print 'Completed: CryoSat Data Read!'
    
    return thickness_cryosat,lat,lon
    
def joinGrids(thickness_icesat,thickness_cryosat,lat,lon):

    year2010 = np.empty((1,lat.shape[0],lat.shape[1]))
    year2010.fill(np.nan)
    
    sit_s = np.append(thickness_icesat,year2010,axis=0)  
    sit_s = np.append(sit_s,thickness_cryosat,axis=0) 
    
    print "Completed: Calculate Climatology!"    
    return sit_s
    
def netcdfSatelliteJ(lats,lons,var,directory):
    name = 'satellite_regrid_March_20042015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Satellite data for ICESat-J (2004-2009) and' \
                         'CryoSat-2 (2011-2015) that have been regridded' \
                         'on 180x180 EASE2.0 100km grids. Year 2010 is' \
                         'available in the data but has been filled with' \
                         'nan values to be inclusive. Shape is therefore' \
                         '[12,180,180]'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('thick','f4',('years','lat','lon'))
    
    ### Metrics
    varns.units = 'meters'
    ncfile.title = 'ICESat-J/CryoSat March SIT'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NASA-J/ESA Products'
    ncfile.references = 'R. Kwok, S. Laxon'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'
    
### Call functions
sit_ij,lat_ij,lon_ij = readICESatJ()
sit_c,lat_c,lon_c = readCryoSat()
#sit_s = joinGrids(sit_ij,sit_c,lat_ij,lon_ij)
#netcdfSatelliteJ(lat_ij,lon_ij,sit_s,directorynew)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
    
def readICESatJ(directory): 
    filename = directory + 'icesatG_regrid_March_20032008.nc'
    data = Dataset(filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    thickness_icesat = data.variables['sit'][:]
    data.close()
    
    print 'Completed: ICESat-G Data Read!'
    
    return thickness_icesat,lat,lon
    
def joinGrids(thickness_icesat,thickness_cryosat,lat,lon):

    year20092010 = np.empty((2,lat.shape[0],lat.shape[1]))
    year20092010.fill(np.nan)
    
    sit_s = np.append(thickness_icesat,year20092010,axis=0)  
    sit_s = np.append(sit_s,thickness_cryosat,axis=0) 
    
    print "Completed: Calculate Climatology!"    
    return sit_s
    
def netcdfSatelliteG(lats,lons,var,directory):
    name = 'satelliteG_regrid_March_20032015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Satellite data for ICESat-G (2003-2008) and' \
                         'CryoSat-2 (2011-2015) that have been regridded' \
                         'on 180x180 EASE2.0 100km grids. Years 2009-2010 are' \
                         'available in the data but have been filled with' \
                         'nan values to be inclusive. Shape is therefore' \
                         '[13,180,180]'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('thick','f4',('years','lat','lon'))
    
    ### Metrics
    varns.units = 'meters'
    ncfile.title = 'ICESat_G/CryoSat March SIT'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NASA-G/ESA Products'
    ncfile.references = 'Donghui Yi, H. Zwally, S. Laxon'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'

### Call Functions
sit_ig,lat_ig,lon_ig = readICESatJ(directorynew)
sit_sg = joinGrids(sit_ig,sit_c,lat_ig,lon_ig)
#netcdfSatelliteG(lat_ig,lon_ig,sit_sg,directorynew)
