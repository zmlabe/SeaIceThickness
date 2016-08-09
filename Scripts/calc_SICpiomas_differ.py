"""
Script compares data between sea ice concentration data from NSIDC versus
PIOMAS sea ice concentration. Monthly time scales.
 
Author : Zachary Labe
Date : 13 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime

### Define directories
directorydata1 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/SeaIceConcentration/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/sic/'  
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1993
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- NSIDC & PIOMAS Comparisons (%s) ---' '\n' % titletime 

def satReader(directory,month):
    """
    Reads NSIDC sea ice concentration data (regridded), 1993-2015
    """
    
    ### Enter filename
    filename = 'nsidc_regrid_sic_19932015.nc'    
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    sic_s = data.variables['sic'][:,month,:,:]
    data.close()
    
    print 'Completed: Satellite data read (March)!'
    
    return lats,lons,sic_s

def piomasReader(directory,month,years):
    """
    Reads piomas data for sea ice concentration over 1979-2015
    """
    
    ### Enter filename
    filename = 'piomas_regrid_sic_19792015.nc'   
    
    ### Month/Years extracted
    dateyr = now.year  
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    ### Retrieve data
    data = Dataset(directory + filename)
    latp = data.variables['lat'][:]
    lonp = data.variables['lon'][:]
    sic_p = data.variables['sic'][:,month,:,:]
    data.close()
    
    sic_p[np.where(sic_p == 1.0)] = np.nan
    
    print 'Completed: PIOMAS data read (%s)!' % datemo
    
    return latp,lonp,datemo,sic_p
    
def diff(sic_s,sic_p):
    """
    Subtracts satellite-piomas
    """
    sic_diff = sic_s[:,:,:] - sic_p[-23:,:,:]
    
    print 'Completed: Difference calculation!'
    
    return sic_diff
    
def plotSubplot(lats,lons,datemo,year,years,sic_s,sic_p,grid,directory):
    """
    Plots subplot 1x2 for sea ice concentration for ortho or polar grids
    """
    
    dateyr = years[year]
    
    def colormapSIC():
        cmap = plt.get_cmap('RdPu')
        cmaplist = [cmap(i) for i in xrange(0,cmap.N-20)]
        cms_sic = c.ListedColormap(cmaplist)
    
        return cms_sic
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    
    ### Set data limits
    lim = np.arange(0.2,1.05,0.05) 
    sic_s[np.where(sic_s < 0.15)] = 0.15
    sic_p[np.where(sic_p < 0.15)] = 0.15
    info = 'Sea Ice Concentration'
    
    ### Plot filled contours 
    ax_s = plt.subplot(121)   
    ### Set grid style
    if grid == 'ortho':
        m = Basemap(projection='ortho',lon_0=-90,
                    lat_0=70,resolution='l',round=True)
    elif grid == 'polar':
        m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
                    resolution='l',round =True)
            
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.4,color='k',fontsize=4)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.4,color='k',fontsize=4)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    cs_s = m.contourf(lons[:,:],lats[:,:],sic_s[year,:,:],
                    lim,latlon=True,extend='max')
    cs_s.set_cmap(colormapSIC())
###########################################################################
###########################################################################    
    ax_p = plt.subplot(122)
    ### Set grid style
    if grid == 'otho':
        m = Basemap(projection='ortho',lon_0=-90,
                    lat_0=70,resolution='l',round=True)
    elif grid == 'polar':
        m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
                    resolution='l',round =True)
            
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.4,color='k',fontsize=4)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.4,color='k',fontsize=4)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    cs_p = m.contourf(lons[:,:],lats[:,:],sic_p[year,:,:],
                lim,latlon=True,extend='max',drawedges=True)
    cs_p.set_cmap(colormapSIC())
    
    cbar_ax = fig.add_axes([0.3109,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs_s,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac='auto')

    cbar.set_label(r'Thickness (meters)')
    cbar.set_ticks(np.arange(0.2,1.1,.1))
    cbar.set_ticklabels(map(str,np.arange(0.2,1.1,.1)))       
    
    fig.suptitle(r'\textbf{%s %s %s}' % (info,datemo,dateyr),
                 fontsize=16)
    plt.annotate(r'\textbf{PIOMAS}',xy=(1.2,21),xycoords='axes fraction',
         fontsize=22)
    plt.annotate(r'\textbf{NSIDC}',xy=(-0.5,21),
                 xycoords='axes fraction',fontsize=26)
#    fig.subplots_adjust(top=0.85)
    fig.subplots_adjust(wspace=0.3)
    
    print 'Completed: Subplot finished!'
                
    plt.savefig(directory +'SICpiomas_%s_%s.png' % (datemo,dateyr),
                dpi=500)
                
def plotDifference(lats,lons,years,sic_diff,directory):
    """
    Makes large subplot of difference in all the SIC data between
    PIOMAS and the satellites [satellite - PIOMAS]
    """
        
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    
    for i in xrange(len(sic_diff)):
        ax = plt.subplot(4,6,i+1)    
        m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,
                    round=True,resolution='l')
        m.drawmapboundary(fill_color = 'white')
        m.drawlsmask(land_color='darkgrey',ocean_color='snow')
        parallels = np.arange(50,90,10)
        meridians = np.arange(-180,180,30)
        m.drawparallels(parallels,labels=[False,False,False,False],
                        linewidth=0.25)
        m.drawmeridians(meridians,labels=[False,False,False,False],
                        linewidth=0.25)

        sic_diff[np.where(sic_diff >= 1)] = 1
        sic_diff[np.where(sic_diff <= -1)] = -1
        
        cs = m.contourf(lons,lats,sic_diff[i,:,:],np.arange(-1,1.01,0.01),
                        latlon=True,extend='both')
        cs.set_cmap('seismic_r')
        ax.text(0.89,0.95,r'\textbf{%s}' % (years[i]),size='5',
                horizontalalignment='center',backgroundcolor='w',
                verticalalignment='center',bbox=dict(facecolor='w',
                edgecolor='k',alpha=0.9),transform=ax.transAxes)
                
    cbar_ax = fig.add_axes([0.30,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac=0.07)

    cbar.set_label(r'SIC ($\times$100\%)')
    cbar.set_ticks(np.arange(-1,1.1,.5))
    cbar.set_ticklabels(map(str,np.arange(-1,1.1,.5)))   
    
    fig.suptitle(r'\textbf{SIC Difference [Satellite - PIOMAS]}',fontsize=16)
    plt.tight_layout() 
    fig.subplots_adjust(bottom=0.14)
    fig.subplots_adjust(top=0.92)
    fig.subplots_adjust(wspace=-.77)         
                
    plt.savefig(directory +'SICpiomas_diff.png',
                dpi=500)
    print 'Completed: Difference subplot finished!'
        
### Call functions
lats,lons,sic_s = satReader(directorydata2,2)  
latp,lonp,datemo,sic_p = piomasReader(directorydata1,2,years) 
sic_diff = diff(sic_s,sic_p)  

plotSubplot(lats,lons,datemo,-2,years,sic_s,sic_p,'polar',directoryfigure)
#plotDifference(lats,lons,years,sic_diff,directoryfigure)

print 'Completed: Script done!'