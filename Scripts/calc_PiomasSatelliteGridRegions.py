"""
Script creates scatter plots between March data for SIT using PIOMAS
and ICESat-J and CryoSat-2. Subplots are for each individual year 
over the 2004-2015 time frame. Scripts makes grid cells larger to see how
smooth (or not) the distribution of SIT is over the Arctic.
 
Author : Zachary Labe
Date : 18 August 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import iris as ir
import iris.quickplot as qplt

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
yearssub = np.arange(1986,1994+1,1)
timex = np.arange(0,8,1)
timey = np.arange(0,8,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- PIOMAS/Sat Yearly Scatter Plots & Large Cells (%s) ---' '\n' % titletime 

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
    else:
        sitp = data.variables['thick'][:,:,:]
        
    data.close()
    
    print 'Completed: PIOMAS data read!'
    return lat,lon,sitp
    
def icesatReader(directory):
    
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
      
### Call functions
lat,lon,sitpi = piomasReader(directorydata,'icej',years)
lat,lon,sitpc = piomasReader(directorydata,'cryo',years)
lat,lon,siti,sitc = icesatReader(directorydata)

def transformGrid(var,la,lo,var2,lat,lon,types):
    
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
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
    
    if types == 'icesat':
        typed = 'ICESatRegions/'
        years = ['2004/','2005/','2006/','2007/','2008/','2009/']
        labely = 'sit( ICESat-J )'
        color = 'seagreen'
        edgecolor = 'darkgreen'
    elif types == 'cryosat':
        typed = 'CryoSat2Regions/'
        years = ['2011/','2012/','2013/','2014/','2015/']
        labely = 'sit( CryoSat-2 )'
        color = 'steelblue'
        edgecolor = 'darkblue'
        
    varn_re = np.empty(var.shape)
    for i in xrange(var.shape[0]):
        
        directoryfigure2 = '/home/zlabe/Desktop/LargeGridCells/' + typed + years[i]
        for j in xrange(0,var.shape[1]-la,la):
            for k in xrange(0,var.shape[2]-lo,lo):
                averaging = np.nanmean(var[i,j:j+la,k:k+lo])
                averaging2 = np.nanmean(var2[i,j:j+la,k:k+lo])
                varn_re[i,j:j+la,k:k+lo] = averaging
                
                if np.isfinite(averaging):
                    if np.isfinite(averaging2):
                        
                        fig = plt.figure()
                        ax = plt.subplot(111)
                        
                        plt.scatter(var[i,j:j+la,k:k+lo],var2[i,j:j+la,k:k+lo],
                                    color='seagreen',s=35,zorder=3,
                                    edgecolor='darkgreen',linewidth=0.5)
                        plt.plot(timex,timey,linewidth=2,color='k',zorder=1)
                        
                        varx = np.ravel(var[i,j:j+la,k:k+lo])
                        vary = np.ravel(var2[i,j:j+la,k:k+lo])
                        mask = np.isfinite(varx) & np.isfinite(vary)           
                        
                        rss = np.empty((var.shape[0],var.shape[1],var.shape[2]))
                        if any(mask):
                            fit = np.polyfit(varx[mask],vary[mask],1)
                            m = fit[0]
                            b = fit[1]
                            line = m*timex + b                                
                            plt.plot(line,linewidth=1,color='r',zorder=2)
                            
                            slope,intercept,r,p_value,std_err = sts.stats.linregress(varx[mask],
                                                          vary[mask])
                            ax.text(6,1,r'r$^2$= %s' % round(r,2),color='k',fontsize=15)
                            
                            r = round(r,2)
                            if r > 1.:
                                r = 1.
                            elif r < -1:
                                r = -1.
                            rss[i,j:j+la,k:k+lo] = r
                        else:
                            rss[i,j:j+la,k:k+lo] = np.nan
                        
                        plt.xlim([0,7])
                        plt.ylim([0,7])
                        plt.xticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
                        plt.yticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
                        plt.xlabel(r'\textbf{sit( PIOMAS )}',fontsize=11)
                        plt.ylabel(r'\textbf{%s}' % labely,fontsize=11)
                        ax.spines['top'].set_color('none')
                        ax.spines['right'].set_color('none')
                        adjust_spines(ax, ['left', 'bottom'])
                        ax.tick_params(labeltop='off', labelright='off')
                        ax.xaxis.set_ticks_position('bottom')
                        ax.yaxis.set_ticks_position('left')
                        fig.suptitle(r'\textbf{%s}' % years[0][0:4],fontsize=17)
                        
                        a2 = plt.axes([.15, .66, .29, .29], axisbg='w')   
                        m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
                        m.drawmapboundary(fill_color = 'white')
                        m.drawcoastlines(color = 'darkgrey',linewidth=0.2)
                        m.drawlsmask(land_color='darkgrey',ocean_color='snow')
                        parallels = np.arange(50,90,10)
                        meridians = np.arange(-180,180,30)
                        m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.25)
                        m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.25,
                                        fontsize=4)
                                        
                        filled = var[i,j:j+la,k:k+lo].copy()
                        filled.fill(1)
                        cs = m.contourf(lon[j:j+la,k:k+lo],lat[j:j+la,k:k+lo],filled,np.arange(0,8,1),
                                        latlon=True,colors='seagreen')    
                        m.fillcontinents(color='darkgrey')
                        
                        fig.subplots_adjust(bottom=0.15)
#                        plt.savefig(directoryfigure2 + '%s_%s_%s.png' % (types,j,k),dpi=300)
    return varn_re,rss

#sitnewpi,ri = transformGrid(sitpi,7,8,siti,lat,lon,'icesat')
#sitnewc,rc = transformGrid(sitpc,7,8,sitc,lat,lon,'cryosat') 

#r = ri
r = rc

lat = lat
lon = lon

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
ax = plt.subplot(111)

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
#values = np.arange(0,7.1,0.1)  
#sit[np.where(sit >= 7.0)] = 7.0
#info = 'Sea Ice Thickness'
values = np.arange(-1,1.1,.1)

### Plot filled contours    
cs = m.contourf(lon[:,:],lat[:,:],r[0,:,:],
                values,latlon=True,extend='max')
                  
### Set colormap     
cmap = plt.cm.get_cmap('brewer_RdBu_11')                           
cs.set_cmap(cmap)

### Set colorbar
cbar = m.colorbar(cs,drawedges=True,location='right',pad = 0.55)

cbar.set_ticks(np.arange(-1,1.5,0.5))
cbar.set_ticklabels(map(str,np.arange(-1,1.5,0.5)))    
cbar.set_label(r'\textbf{r$^2$}')

fig.suptitle(r'\textbf{ICESat-J}',fontsize=13)    
fig.subplots_adjust(top=0.85)

### Save figure
plt.savefig(directoryfigure +'rsquared.png',dpi=500)