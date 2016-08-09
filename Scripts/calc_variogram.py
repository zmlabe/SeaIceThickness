"""
Testing scripts for calculating variograms with PIOMAS data
 
Author : Zachary Labe
Date : 5 August 2016
Source : pydoc.net/Python/ambhas/0.4.0/ambhas.krige/
"""

### Import Modules
import datetime
import numpy as np
from netCDF4 import Dataset
from scipy.spatial.distance import pdist,squareform
import matplotlib.pyplot as plt

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

print '\n' '--- PIOMAS Variogram (%s) ---' '\n' % titletime 

### Read in data
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
    elif segment == 'all':
        sitp = data.variables['thick'][:,:,:]
        
    data.close()
    
    print 'Completed: PIOMAS data read!'
    return lat,lon,sitp
    
### Use 1 year to get 2d array    
lat,lon,sitp = piomasReader(directorydata,'all',years)
sitp = sitp[-1,:,:]

x = lat
y = lon
z = sitp
D = np.sqrt(lon**2 + lat**2)
g = 0.5*z**2
indx = range(len(z))
C,R = np.meshgrid(indx,indx)
I = R > C

D2 = D*(np.diag(x*np.nan)+1)
lag = 10*np.mean(np.nanmin(D2))