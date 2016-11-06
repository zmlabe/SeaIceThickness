"""
Functions are reader for March sea ice thickness data.
 
Author : Zachary Labe
Date : 11 August 2016
"""

def piomasReader(segment):
    
    import numpy as np
    from netCDF4 import Dataset
    
    yearmin = 1979
    yearmax = 2015
    years = np.arange(yearmin,yearmax+1,1)
    
    directory = '/home/zlabe/Surtsey/seaice_obs/Thk/March/' 
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
    
def icesatReader():
    
    from netCDF4 import Dataset
    
    directory = '/home/zlabe/Surtsey/seaice_obs/Thk/March/' 
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