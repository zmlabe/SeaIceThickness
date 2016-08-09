"""
Create climatology for ICESat and CryoSat-2 data by joining grids
 
Source : ftp://sidads.colorado.edu/pub/projects/SIPN/seaice_thickness/
Author : Zachary Labe
Date : 12 July 2016
"""

### Import Modules
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

### Define constants
year_icesat = np.arange(2004,2010)
year_cs2 = np.arange(2011,2016)
years = np.append(year_icesat,year_cs2)

# Read in Data
def icesat(year):   
    directory = '/home/zlabe/Documents/SeaIceObs/Data/'
    
    thickness_icesat = np.empty((len(year),180,180))
    for i in xrange(len(year)):
        files = 'icesat_seaice_thickness_mean_spring_%s.nc' % year[i]
        filename = directory + files
        data1 = Dataset(filename)
        lat = data1.variables['lat'][:]
        lon = data1.variables['lon'][:]
        thickness_icesat[i,:,:] = data1.variables['thick'][:]
        data1.close()
        
    thickness_icesat[np.where(thickness_icesat == 9999.)] = np.nan
    lon[np.where(lon == -999.)] = np.nan
    lat[np.where(lat == -999.)] = np.nan
    
    print 'Completed: IceSat Data Read!'
        
    return thickness_icesat,lat,lon
    
def cs2():
    directory = '/home/zlabe/Documents/SeaIceObs/Data/'
    
    filename = directory + 'cryosat_March_20112015.nc'
    data2 = Dataset(filename)
    lat2 = data2.variables['lat'][:]
    lon2 = data2.variables['lon'][:]
    thickness_cryosat = data2.variables['newthickness'][:]
    data2.close()
    
    print 'Completed: CryoSat Data Read!'
    
    return thickness_cryosat,lat2,lon2
    
def joinGrids(thickness_icesat,thickness_cs2,lat,lon):
    thickness = np.append(thickness_icesat,thickness_cs2,axis=0)  
    
    print "Completed: Calculate Climatology!"    
    return thickness
    
def netcdf(lats,lons,var):
    directory = '/home/zlabe/Surtsey/seaice_obs/Sat_thick/'
    name = 'cs2icesat_regrid_mar_20042015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'Join ICESat and CS-2 on same grid for years' \
                        ' 2004-2015 during month of March'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('lat',var.shape[1])
    ncfile.createDimension('lon',var.shape[2])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('thick','f4',('years','lat','lon'))
    
    ### Units
    varns.units = 'meters'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print 'Completed: Created netCDF4 File!'
    
### Call functions    
thickness_icesat,lat,lon = icesat(year_icesat)
thickness_cs2,lat2,lon2 = cs2()
thickness = joinGrids(thickness_icesat,thickness_cs2,lat,lon)
netcdf(lat,lon,thickness)

print 'Completed: Script Finished!'