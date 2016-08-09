"""
Script plots oceanic heat flux (PIOMAS) data using Basemap module 
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
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/OceanFlux/'    
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1979
yearmax = 2004
years = np.arange(yearmin,yearmax+1,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'Plot Oceanic Heat Flux - %s' '\n' % titletime 

def readPiomas(directory):
    """
    Reads netCDF4 PIOMAS data
    """

    files = 'piomas_regrid_oflux_19792004.nc'        
    
    data = Dataset(directory + files)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    oflux = data.variables['snow'][:]
    data.close()
    
    ### Convert to new units (seconds --> days)
    oflux = oflux * 60 * 60 * 24

    print 'Completed: Read data!'  
    
    return lats,lons,oflux
    
def plotSITsingle(directory,lats,lons,sic,month,year,var,style,thresh):
    """
    Plots sit data for ortho/npstere basemap grids
    """
    
    def colormapOflux():
        cmap = plt.get_cmap('hot_r')
        cmaplist = [cmap(i) for i in xrange(60,cmap.N-40)]
        cms_sic = c.ListedColormap(cmaplist)
        return cms_sic
    
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
    values = np.arange(0,0.1005,0.005)
    oflux[np.where(oflux <= 0.0001)] = 0.0001
    oflux[np.where(oflux >= 0.1)] = 0.1
    info = 'Oceanic Heat Flux'
    
    ### Plot filled contours    
    cs = m.contourf(lons[:,:],lats[:,:],oflux[year,month,:,:],
                    values,latlon=True,extend='max')
                      
    ### Set colormap                              
    cs.set_cmap(colormapOflux())
    
    ### Set colorbar
    cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)
    cbar.set_ticks(np.arange(0,0.15,0.01))
    cbar.set_ticklabels(map(str,np.arange(0,0.15,0.01)))
    cbar.set_label(r'\textbf{Ocean Heat Flux (meters of ice/day)}')
    
    ### Date and title/adjustments
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
        
    fig.suptitle(r'\textbf{%s %s %s [PIOMAS]}' % (datemo,dateyr,info),
             fontsize=15)
    fig.subplots_adjust(top=0.85)
    
    ### Save figure
    plt.savefig(directory +'%s_%s%s_regrid.png' % (var,datemo[0:3],dateyr),dpi=500)

### Call functions
lats,lons,oflux = readPiomas(directorydata)
plotSITsingle(directoryfigure,lats,lons,oflux,11,-7,'oflux','ortho',0.0)

print 'Completed: Script finished!'