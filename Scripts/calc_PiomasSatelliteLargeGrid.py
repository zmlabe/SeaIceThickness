"""
Scripts creates plots of large grid cells (nxn) for different statistical
variables. Other scripts are available for looking at individual
scatter plots over these larger domains.
 
Author : Zachary Labe
Date : 2 September 2016
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
yearsi = np.arange(2004,2010,1)
yearsc = np.arange(2011,2016,1)
timex = np.arange(0,8,1)
timey = np.arange(0,8,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- PIOMAS/Sat Large Cell Statistics (%s) ---' '\n' % titletime 

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
      
    ms = np.zeros((var.shape[0],var.shape[1],var.shape[2]))
    bs = np.zeros((var.shape[0],var.shape[1],var.shape[2]))
    rs = np.zeros((var.shape[0],var.shape[1],var.shape[2]))
    varn_re = np.empty(var.shape)
    for i in xrange(var.shape[0]):
        for j in xrange(0,var.shape[1]-la,la):
            for k in xrange(0,var.shape[2]-lo,lo):
                averaging = np.nanmean(var[i,j:j+la,k:k+lo])
                averaging2 = np.nanmean(var2[i,j:j+la,k:k+lo])
                varn_re[i,j:j+la,k:k+lo] = averaging
                
                if np.isfinite(averaging):
                    if np.isfinite(averaging2):
                        
                        varx = np.ravel(var[i,j:j+la,k:k+lo])
                        vary = np.ravel(var2[i,j:j+la,k:k+lo])
                        mask = np.isfinite(varx) & np.isfinite(vary)           
                        
                        if any(mask):
                            fit = np.polyfit(varx[mask],vary[mask],1)
                            m = fit[0]
                            b = fit[1]                              
                            
                            slope,intercept,r,p_value,std_err = sts.stats.linregress(varx[mask],
                                                          vary[mask])
                            m = round(m,2)
                            ms[i,j:j+la,k:k+lo] = m
                            b = round(b,2)
                            bs[i,j:j+la,k:k+lo] = b
                            r = round(r,2)
                            rs[i,j:j+la,k:k+lo] = r
                else:
                    ms[i,j:j+la,k:k+lo] = np.nan
                    bs[i,j:j+la,k:k+lo] = np.nan
                    rs[i,j:j+la,k:k+lo] = np.nan
                        
    return varn_re,ms,bs,rs

sitnewpi,mi,bi,ri = transformGrid(sitpi,4,5,siti,lat,lon,'icesat')
sitnewc,mc,bc,rc = transformGrid(sitpc,4,5,sitc,lat,lon,'cryosat') 

mm = mi
#mm= bc

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

### Define figure
fig = plt.figure()
for i in xrange(mm.shape[0]):
    ax = plt.subplot(2,3,i+1)
    
    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.5)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0.5,color='k',fontsize=9)
#    m.drawmeridians(meridians,labels=[True,True,False,False],
#                    linewidth=0.5,color='k',fontsize=9)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(-2,2.1,0.25)  
#    values = np.arange(-5,5.1,0.25)
#    values = np.arange(-1,1.1,0.2)
    
    mm[np.where(mm == 0.)] = np.nan
    
    ### Plot filled contours    
    cs = m.contourf(lon[:,:],lat[:,:],mm[i,:,:],
                    values,latlon=True,extend='both')
                      
    ### Set colormap     
#    cmap = plt.cm.get_cmap('brewer_RdBu_11')      
    cmap = plt.cm.get_cmap('bwr')                      
    cs.set_cmap(cmap)
    
    ### Set colorbar
    ax.text(0.92,0.97,r'%s' % (yearsi[i]),size='13',
            horizontalalignment='center',
            verticalalignment='center',transform=ax.transAxes)
            
cbar_ax = fig.add_axes([0.31,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'slopes (m)')
cbar.set_ticks(np.arange(-2,2.1,0.5))
cbar.set_ticklabels(map(str,np.arange(-2,2.1,0.5))) 
#cbar.set_ticks(np.arange(-5,5.1,2))
#cbar.set_ticklabels(map(str,np.arange(-5,5.1,2))) 
#cbar.set_ticks(np.arange(-1,1.1,0.5))
#cbar.set_ticklabels(map(str,np.arange(-1,1.1,0.5))) 
cbar_ax.tick_params(axis=u'both', which=u'both',length=0)

fig.subplots_adjust(bottom=0.15)
plt.text(0.1,28,r'\textbf{PIOMAS - ICESat-J}',fontsize=15)

### Save figure
plt.savefig(directoryfigure +'icesat_slopes.png',dpi=500)