"""
Script plots snow cover (PIOMAS) data using Basemap module 
between ortho or polar stereographic grids
 
Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
Author : Zachary Labe
Date : 11 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/SnowCover/'    
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

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

print '\n' 'Plot Snow Cover - %s' '\n' % titletime 

def readPiomas(directory):
    """
    Reads netCDF4 PIOMAS data
    """

    files = 'piomas_regrid_snow_19792015.nc'        
    
    data = Dataset(directory + files)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    sit = data.variables['snow'][:]
    data.close()

    print 'Completed: Read data!'  
    
    return lats,lons,sit
    
def plotSITsingle(directory,lats,lons,snow,month,year,var,style,thresh):
    """
    Plots sit data for ortho/npstere basemap grids
    """
    
    def colormapSnow():
        cmap1 = plt.get_cmap('viridis')
        cms_snow = cmap1 
        return cms_snow
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    ax = plt.subplot(111)
    
    if style == 'ortho':
        m = Basemap(projection='ortho',lon_0=-90,
                    lat_0=70,resolution='l',round=True)
    elif style == 'polar':
        m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
        
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.5,color='k',fontsize=9)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.5,color='k',fontsize=9)
    m.drawlsmask(land_color='peru',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(0,0.31,0.01)
#    snow[np.where(snow >= 0.3)] = 0.3
    info = 'Snow Cover'
    
    ### Plot filled contours    
    cs = m.contourf(lons[:,:],lats[:,:],snow[year,month,:,:],
                    values,latlon=True,extend='max')
                      
    ### Set colormap                              
    cs.set_cmap(colormapSnow())
    
    ### Set colorbar
    cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)
    
    cbar.set_ticks(np.arange(0.0,0.31,0.05))
    cbar.set_ticklabels(map(str,np.arange(0.0,0.31,0.05)))
    cbar.set_label(r'\textbf{Snow Depth (meters)}')
    
    ### Date and title/adjustments
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
        
    fig.suptitle(r'\textbf{%s %s %s [PIOMAS]}' % (datemo,dateyr,info),
             fontsize=15)
    fig.subplots_adjust(top=0.85)
    
    ### Save figure
    plt.savefig(directory +'%s_%s%s_regrid.png' % (var,datemo[0:3],dateyr),dpi=500)

### Call functions
lats,lons,snow = readPiomas(directorydata)
plotSITsingle(directoryfigure,lats,lons,snow,2,0,'snow','ortho',0.0)

print 'Completed: Script finished!'