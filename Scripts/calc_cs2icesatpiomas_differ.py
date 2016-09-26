"""
Script compares data between altimetric satellites (CS-2 & ICESat) with
data from PIOMAS. Functions for zooming on particular grid cell locations
in addition to the entire Arctic basin. 
 
Author : Zachary Labe
Date : 12 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime

### Define directories
directorydata1 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/Sat_thick/'  
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 2004
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
years = np.setdiff1d(years,[2010])      ### no satellite data in 2010
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- CS2/ICESat & PIOMAS Comparisons (%s) ---' '\n' % titletime 

def satReader(directory):
    """
    Reads satellite data for joint CS-2/ICESat gridded (2004-2015)
    """
    
    ### Enter filename
    filename = 'cs2icesat_regrid_mar_20042015.nc'    
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    thk_s = data.variables['thick'][:]
    data.close()
    
    print 'Completed: Satellite data read (March)!'
    
    return lats,lons,thk_s

def piomasReader(directory,month,years):
    """
    Reads piomas data for sea ice thickness over 1979-2015
    """
    
    ### Enter filename
    filename = 'piomas_regrid_sit_19792015.nc'   
    
    ### Month/Years extracted
    dateyr = now.year  
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    yearsp = np.arange(1979,2016)
    yearmin = years.min()
    yearmax = years.max()
    yearnone = 2010
    yearslice = np.where((yearsp <= yearmax) & (yearsp >= yearmin) & \
                (yearsp != yearnone))[0]
    
    ### Retrieve data
    data = Dataset(directory + filename)
    latp = data.variables['lat'][:]
    lonp = data.variables['lon'][:]
    thk_p = data.variables['newthickness'][yearslice,month,:,:]
    data.close()
    
    print 'Completed: PIOMAS data read (%s)!' % datemo
    return latp,lonp,thk_p
    
def diff(thk_s,thk_p):
    """
    Subtracts satellite-piomas
    """
    thk_diff = thk_s - thk_p
    
    print 'Completed: Difference calculation!'
    return thk_diff
    
def plotSubplot(lats,lons,year,years,thk_s,thk_p,grid,directory):
    """
    Plots subplot 1x2 for thickness for ortho or polar grids
    """
    
    def colormapSIT():
        cmap1 = plt.get_cmap('BuPu')
        cmap2 = plt.get_cmap('RdPu_r')
        cmap3 = plt.get_cmap('gist_heat_r')
        cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
        cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
        cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
        cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
        return cms_sit
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    
    ### Set data limits
    lim = np.arange(0,7.1,0.1)  
    thk_s[np.where(thk_s >= 7.0)] = 7.0 
    thk_p[np.where(thk_p >= 7.0)] = 7.0
    info = 'Sea Ice Thickness'
    
    ### Plot filled contours 
    ax_s = plt.subplot(121)   
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
                    linewidth=0.5,color='k',fontsize=4)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.5,color='k',fontsize=4)
    m.drawlsmask(land_color='peru',ocean_color='mintcream')
    
    cs_s = m.contourf(lons[:,:],lats[:,:],thk_s[year,:,:],
                    lim,latlon=True,extend='max')
    cs_s.set_cmap(colormapSIT())
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
                    linewidth=0.5,color='k',fontsize=4)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.5,color='k',fontsize=4)
    m.drawlsmask(land_color='peru',ocean_color='mintcream')
    cs_p = m.contourf(lons[:,:],lats[:,:],thk_p[year,:,:],
                lim,latlon=True,extend='max')
    cs_p.set_cmap(colormapSIT())
    
    dateyr = years[year]    
    
    if dateyr > 2009:
        satellite = 'CS-2'
    elif dateyr < 2011:
        satellite = 'ICESat'
    
    cbar_ax = fig.add_axes([0.3109,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs_s,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac='auto')

    cbar.set_label(r'Thickness (meters)')
    cbar.set_ticks(np.arange(0,7.1,1))
    cbar.set_ticklabels(map(str,np.arange(0,7.1,1)))       
    
    fig.suptitle(r'\textbf{Sea Ice Thickness March %s}' % dateyr,
                 fontsize=16)
    plt.annotate(r'\textbf{PIOMAS}',xy=(1.2,21),xycoords='axes fraction',
         fontsize=22)
    plt.annotate(r'\textbf{%s}' % satellite,xy=(-0.5,21),
                 xycoords='axes fraction',fontsize=26)
#    fig.subplots_adjust(top=0.85)
    fig.subplots_adjust(wspace=0.3)
    
    print 'Completed: Subplot finished!'
                
    plt.savefig(directory +'%spiomas_mar_%s.png' % (satellite,dateyr),
                dpi=500)
                
def plotDifference(lats,lons,years,thk_diff,directory):
    """
    Makes large subplot of difference in all the SIT data between
    PIOMAS and the satellites [satellite - PIOMAS]
    """
        
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    
    for i in xrange(len(thk_diff)):
        ax = plt.subplot(3,4,i+1)    
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

        thk_diff[np.where(thk_diff >= 3.)] = 3
        thk_diff[np.where(thk_diff <= -3.)] = -3
        
        cs = m.contourf(lons,lats,thk_diff[i,:,:],np.arange(-3,3.1,0.1),
                        latlon=True,extend='both')
        cs.set_cmap('seismic_r')
        ax.text(0.89,0.95,r'\textbf{%s}' % (years[i]),size='8',
                horizontalalignment='center',backgroundcolor='w',
                verticalalignment='center',bbox=dict(facecolor='w',
                edgecolor='k',alpha=0.9),transform=ax.transAxes)
                
    cbar_ax = fig.add_axes([0.30,0.1,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac=0.07)

    cbar.set_label(r'Thickness (m)')
    cbar.set_ticks(np.arange(-3,4,1))
    cbar.set_ticklabels(map(str,np.arange(-3,4,1)))   
    
    fig.suptitle(r'\textbf{SIT Difference [Satellite - PIOMAS]}',fontsize=16)
    plt.tight_layout() 
    fig.subplots_adjust(bottom=0.14)
    fig.subplots_adjust(top=0.92)
    fig.subplots_adjust(wspace=-.56)
    fig.subplots_adjust(hspace=-0.0)         
                
    plt.savefig(directory +'cs2piomas_diff.png',
                dpi=500)
    print 'Completed: Difference subplot finished!'
        
        
### Call functions
lats,lons,thk_s = satReader(directorydata2)  
latp,lonp,thk_p = piomasReader(directorydata1,2,years) 
thk_diff = diff(thk_s,thk_p)  

plotSubplot(lats,lons,-10,years,thk_s,thk_p,'polar',directoryfigure)
plotDifference(lats,lons,years,thk_diff,directoryfigure)
#
#
#fig = plt.figure()
#
#def colormapSIT():
#    cmap1 = plt.get_cmap('BuPu')
#    cmap2 = plt.get_cmap('RdPu_r')
#    cmap3 = plt.get_cmap('gist_heat_r')
#    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
#    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
#    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
#    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
#    return cms_sit


#a = np.nanmean(thk_s,axis=2)
#d = thk_s - thk_p
#d = np.transpose(np.nanmean(d,axis=2))
#d[np.where(d >= 3)] = 3.
#d[np.where(d <= -3)]= 3.
#
#a = np.transpose(a)
#
#c1 = plt.contourf(d,np.arange(-3,3.1,0.1))
#
#labelsx = map(str,years)
#plt.xticks(np.arange(len(years)),labelsx)
##plt.yticks([])
##c1.set_cmap(colormapSIT())
#c1.set_cmap('seismic_r')
#cb = plt.colorbar(c1,extend='both')
#plt.axvline(6,color='k',linewidth=3,linestyle='-')
#cb.set_ticks(np.arange(-3,4,1))
#cb.set_ticklabels(map(str,np.arange(-3,4,1)))
#cb.set_label(r'Thickness (meters)')
#plt.ylim([50,120])
#fig.suptitle(r'\textbf{Zonal Mean Difference [Satellite - PIOMAS]}',
#             fontsize = 15)
#
#plt.savefig(directoryfigure + 'testhov.png',dpi=300)

print 'Completed: Script done!'