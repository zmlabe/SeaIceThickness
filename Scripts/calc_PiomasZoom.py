"""
Script zooms in on small grid cell regions for PIOMAS sea ice thickness
using function with defined latitude x longitude regions. Script creates 
gridded plot and lines plots for comparison.
 
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
import math

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/Sat_thick/'  
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

print '\n' '--- PIOMAS Zoomed Regions (%s) ---' '\n' % titletime 

def piomasReader(directory,month,latmin,latmax,lonmin,lonmax):
    """
    Reads piomas data for sea ice thickness over 1979-2015
    """
    
    ### Enter filename
    filename = 'piomas_regrid_sit_19792015.nc'   
    
    ### Month/Years extracted
    dateyr = now.year  
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    thkn = data.variables['newthickness'][-11:,month,:,:]
    data.close()
    
    ### Calculate lat/lon region
    xmask = (lat > latmin) & (lat < latmax)
    ymask = (lon > lonmin) & (lon < lonmax)
    
    mask = xmask[:] & ymask[:]
    latvals = np.where(mask == True)[0]
    lonvals = np.where(mask == True)[1]
    latvals = np.unique(latvals)
    lonvals = np.unique(lonvals)
    
    thk = thkn[:,latvals,:]
    thk = thk[:,:,lonvals]
    
    lat = lat[latvals,:]
    lat = lat[:,lonvals]
    lon = lon[latvals,:]
    lon = lon[:,lonvals]

    grid = '---> [[%s to %s N, %s to %s E]]' % (latmin,latmax,lonmin,lonmax)
    print 'Completed: PIOMAS data read (%s)!' % datemo, grid
    
    return lat,lon,thk
    
def satReader(directory,month,latmin,latmax,lonmin,lonmax):
    """
    Reads satellite data for joint CS-2/ICESat gridded (2004-2015)
    """
    
    ### Enter filename
    filename = 'cs2icesat_regrid_mar_20042015.nc'  
    
    ### Month/Years extracted
    dateyr = now.year  
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    ### Retrieve data
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    thkn = data.variables['thick'][:]
    data.close()
    
    ### Calculate lat/lon region
    xmask = (lat > latmin) & (lat < latmax)
    ymask = (lon > lonmin) & (lon < lonmax)
    
    mask = xmask[:] & ymask[:]
    latvals = np.where(mask == True)[0]
    lonvals = np.where(mask == True)[1]
    latvals = np.unique(latvals)
    lonvals = np.unique(lonvals)
    
    thk = thkn[:,latvals,:]
    thk = thk[:,:,lonvals]
    
    lat = lat[latvals,:]
    lat = lat[:,lonvals]
    lon = lon[latvals,:]
    lon = lon[:,lonvals]

    grid = '---> [[%s to %s N, %s to %s E]]' % (latmin,latmax,lonmin,lonmax)
    print 'Completed: Satellite data read (%s)!' % datemo, grid
    
    return lat,lon,thk
    
def diff(thk_s,thk_p):
    """
    Subtracts satellite-piomas
    """
    thk_diff = thk_s[:,:,:] - thk_p[:,:,:]
    
    print 'Completed: Difference calculation!'
    return thk_diff

def plotSIT(lat,lon,year,years,month,thk,latmin,latmax,lonmin,lonmax):
    
    latmin = latmin - 0
    latmax = latmax + 0
    lonmin = lonmin - 20
    lonmax = lonmax + 20
    
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
    ax = plt.subplot(121)
    
    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
        
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.4,color='k',fontsize=5)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.4,color='k',fontsize=5)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(0,7.1,0.1)  
    thk[np.where(thk >= 7.0)] = 7.0
    info = 'Sea Ice Thickness'
    
    ### Plot filled contours    
    cs = m.contourf(lon,lat,thk[year,:,:],
                    values,latlon=True,extend='max')
                                  
    ### Set colormap                              
    cs.set_cmap(colormapSIT())
    
    def polar_stere(lon_w, lon_e, lat_s, lat_n, **kwargs):
        '''Returns a Basemap object (NPS/SPS) focused in a region.
        
        lon_w, lon_e, lat_s, lat_n -- Graphic limits in geographical coordinates.
                                      W and S directions are negative.
        **kwargs -- Aditional arguments for Basemap object.
        
        '''
        lon_0 = lon_w + (lon_e - lon_w) / 2.
        ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
        lat_0 = math.copysign(90., ref)
        proj = 'npstere' if lat_0 > 0 else 'spstere'
        prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                              boundinglat=0, resolution='c')
        #prj = pyproj.Proj(proj='stere', lon_0=lon_0, lat_0=lat_0)
        lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
        lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
        x, y = prj(lons, lats)
        ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
        ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
        return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                           llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                           urcrnrlon=ur_lon, urcrnrlat=ur_lat, round=False)      
    
    ax2 = plt.subplot(122)
    m = polar_stere(lonmin,lonmax,latmin,latmax, resolution='l')
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.8)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.4)
    m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.4)
    
    ### Plot filled contours    
    cs = m.contourf(lon,lat,thk[year,:,:],
                    values,latlon=True,extend='max')
    
    ### Set colormap                              
    cs.set_cmap(colormapSIT())
    
    ### Set colorbar
    cbar_ax = fig.add_axes([0.29,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac='auto')

    cbar.set_label(r'Thickness (meters)')
    cbar.set_ticks(np.arange(0,7.1,1))
    cbar.set_ticklabels(map(str,np.arange(0,7.1,1)))  
    
    ### Date and title/adjustments
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
        
    fig.suptitle(r'\textbf{%s %s %s [Satellite]}' % (datemo,dateyr,info),
             fontsize=15)
#    fig.subplots_adjust(top=0.85)
             
    locationgrid = '[%s to %s N, %s to %s E]' % (latmin,latmax,lonmin,lonmax)
    plt.annotate(locationgrid,xy=(0.1,24),xycoords='axes fraction',
         fontsize=13)        
    
    plt.savefig(directoryfigure + 'test5a.png',dpi=300)
    
    print 'Completed: Figure finished!'
    
def plotDifference(lat,lon,year,years,month,thk_diff,latmin,latmax,lonmin,lonmax):
    
    latmin = latmin - 0
    latmax = latmax + 0
    lonmin = lonmin - 20
    lonmax = lonmax + 20

    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define figure
    fig = plt.figure()
    ax = plt.subplot(121)
    
    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
        
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0.4,color='k',fontsize=5)
    m.drawmeridians(meridians,labels=[True,True,False,False],
                    linewidth=0.4,color='k',fontsize=5)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(-3,3.1,0.1)  
    thk_diff[np.where(thk_diff >= 3.)] = 3
    thk_diff[np.where(thk_diff <= -3.)] = -3
    info = 'SIT'
    
    ### Plot filled contours    
    cs = m.contourf(lon,lat,thk_diff[year,:,:],
                    values,latlon=True,extend='both')
                                  
    ### Set colormap                              
    cs.set_cmap('seismic_r')
    
    def polar_stere(lon_w, lon_e, lat_s, lat_n, **kwargs):
        '''Returns a Basemap object (NPS/SPS) focused in a region.
        
        lon_w, lon_e, lat_s, lat_n -- Graphic limits in geographical coordinates.
                                      W and S directions are negative.
        **kwargs -- Aditional arguments for Basemap object.
        
        '''
        lon_0 = lon_w + (lon_e - lon_w) / 2.
        ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
        lat_0 = math.copysign(90., ref)
        proj = 'npstere' if lat_0 > 0 else 'spstere'
        prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                              boundinglat=0, resolution='c')
        #prj = pyproj.Proj(proj='stere', lon_0=lon_0, lat_0=lat_0)
        lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
        lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
        x, y = prj(lons, lats)
        ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
        ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
        return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                           llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                           urcrnrlon=ur_lon, urcrnrlat=ur_lat, round=False)      
    
    ax2 = plt.subplot(122)
    m = polar_stere(lonmin,lonmax,latmin,latmax, resolution='l')
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.8)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.4)
    m.drawmeridians(meridians,labels=[False,False,False,False],linewidth=0.4)
    
    ### Plot filled contours    
    cs = m.contourf(lon,lat,thk_diff[year,:,:],
                    values,latlon=True,extend='both')
    
    ### Set colormap                              
    cs.set_cmap('seismic_r')
    
    ### Set colorbar
    cbar_ax = fig.add_axes([0.29,0.15,0.4,0.03])                
    cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                        extend='both',extendfrac='auto')

    cbar.set_label(r'Thickness (meters)')
    cbar.set_ticks(np.arange(-3,4,1))
    cbar.set_ticklabels(map(str,np.arange(-3,4,1)))  
    
    ### Date and title/adjustments
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
        
    fig.suptitle(r'\textbf{%s %s %s Difference [Satellite - PIOMAS]}' % (datemo,dateyr,info),
             fontsize=15)
#    fig.subplots_adjust(top=0.85)
             
    locationgrid = '[%s to %s N, %s to %s E]' % (latmin,latmax,lonmin,lonmax)
    plt.annotate(locationgrid,xy=(0.1,24),xycoords='axes fraction',
         fontsize=13)        
    
    plt.savefig(directoryfigure + 'test5_diff.png',dpi=300)
    
    print 'Completed: Figure finished!'

def meanRegion(thk_s,thk_p,thk_diff):
    """
    Computes average thickness/differences over defined region
    from previous functions
    """
    meanp = np.nanmean(np.nanmean(thk_p,axis=1),axis=1)
    means = np.nanmean(np.nanmean(thk_s,axis=1),axis=1)
    
    print '\n    --- [[%s to %s N, %s to %s E]] ---' % (latmin,latmax,lonmin,lonmax)
    print 'Average Thickness (Satellite) == %s meters' % np.nanmean(means)
    print 'Average Thickness (PIOMAS)    == %s meters' % np.nanmean(meanp)
    print 'Average Difference            == %s meters' % (np.nanmean(means)-np.nanmean(meanp))
    
    yearmin = 2004
    yearmax = 2015
    years = np.arange(yearmin,yearmax+1,1)
    years = np.setdiff1d(years,[2010])      ### no satellite data in 2010
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    plt.plot(meanp,color='darkred',linewidth=2,linestyle='-',
             label=r'PIOMAS')
    plt.plot(means,color='forestgreen',linewidth=2,linestyle='-',
             label=r'Satellite')
    plt.axvline(6,color='k',linewidth=3,linestyle='-')
    
    labelsy = map(str,np.arange(0,6,1))
    labelsx = map(str,years)
    plt.xticks(np.arange(len(years)),labelsx)
    plt.yticks(np.arange(0,6,1),labelsy)
    plt.ylabel(r'\textbf{Thickness (meters)}',fontsize=13)
    
    ### Adjust axes in time series plots 
    def adjust_spines(ax, spines):
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward', 10))
            else:
                spine.set_color('none')  
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            ax.yaxis.set_ticks([])
    
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            ax.xaxis.set_ticks([]) 
    
    ### Adjust axes spines
    adjust_spines(ax, ['left', 'bottom'])
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    plt.grid(color='b',zorder=1,alpha=0.3)
    
    plt.legend(shadow=False,fontsize=11,loc='upper right',
                       fancybox=True)
    
    plt.text(2,-0.8,r'\textbf{ICESat}',fontsize=13)
    plt.text(7.3,-0.8,r'\textbf{PIOMAS}',fontsize=13)
    
    fig.suptitle(r'\textbf{SIT Difference [Satellite - PIOMAS]}',fontsize=16)
    plt.savefig(directoryfigure + 'test5_difftseries.png',dpi=300)
    
### Call functions
#test 3
latmin = 70
latmax = 89.5
lonmin = 20
lonmax = 90

#test 4
#latmin = 70
#latmax = 89.5
#lonmin = 90
#lonmax = 150

lat,lon,thk_p = piomasReader(directorydata,9,latmin,latmax,lonmin,lonmax)
lat,lon,thk_s = satReader(directorydata2,9,latmin,latmax,lonmin,lonmax) 
#plotSIT(lat,lon,-1,years,9,thk_s,latmin,latmax,lonmin,lonmax)
thk_diff = diff(thk_s,thk_p) 
plotDifference(lat,lon,-1,years,9,thk_diff,latmin,latmax,lonmin,lonmax)
meanRegion(thk_s,thk_p,thk_diff)

print '\nCompleted: Script done!'