"""
Script creates scatter plots between March data for SIT using PIOMAS
and ICESat-J and the submarine data
 
Author : Zachary Labe
Date : 27 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime
import scipy.stats as sts
import matplotlib.mlab as mlab
from matplotlib.ticker import NullFormatter
import statsmodels.api as sm

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
#directoryfigure = '/home/zlabe/Desktop/'
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/' 

yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
yearssub = np.arange(1986,1994+1,1)
       
### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' '--- PIOMAS/Satellites Scatter Plots (%s) ---' '\n' % titletime 

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
    
def subReader(directory):
    
    filename = 'sub_regrid_March_19861994.nc'    
    
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    sitb = data.variables['sit'][:]
    meansitb = data.variables['meansit'][:,:]
    data.close()
    
    print 'Completed: Submarine data read!'
    return lat,lon,sitb,meansitb
    
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
lat,lon,sitp = piomasReader(directorydata,'sub',years)
lat,lon,sitp2 = piomasReader(directorydata,'icej',years)
lat,lon,sitp3 = piomasReader(directorydata,'cryo',years)
lat,lon,sitb,meansitb = subReader(directorydata)
lat,lon,siti,sitc = icesatReader(directorydata)

### Calculate trends
varx = np.ravel(sitp)
vary = np.ravel(sitb)
timex = np.arange(0,8,1)
timey = np.arange(0,8,1)

mask = ~np.isnan(varx) & ~np.isnan(vary)
slope, intercept, r_value, p_value, std_err = sts.linregress(varx[mask], vary[mask])
linep = slope*timex + intercept

#mask = sitb.copy()
#mask[np.where(np.isnan(mask))] = 0.
#mask[np.where(mask != 0.)] = 1.

fig = plt.figure()
ax = plt.subplot(111)

year = 0
month = 2

dateyr = yearssub[year]    
datemo = datetime.date(dateyr,month+1,1).strftime('%B')

### Call parameters
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

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
plt.plot(linep,color='r',zorder=7)

plt.scatter(sitp[0,:,:],sitb[0,:,:],label='1986',color='olive',zorder=2,s=5)
plt.scatter(sitp[1,:,:],sitb[1,:,:],label='1987',color='rosybrown',zorder=3,s=5)
plt.scatter(sitp[5,:,:],sitb[5,:,:],label='1991',color='sandybrown',zorder=4,s=5)
plt.scatter(sitp[-2,:,:],sitb[-2,:,:],label='1993',color='steelblue',zorder=5,s=5)
plt.scatter(sitp[-1,:,:],sitb[-1,:,:],label='1994',color='seagreen',zorder=6,s=5)

### Add legend
plt.legend(shadow=False,fontsize=9,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.93,0.14),
                        frameon=False)
                        
plt.grid(color='darkgrey',linewidth=0.4)

### Add limits
plt.xlim([0,7])
plt.ylim([0,7])
plt.xticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
plt.yticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))

### Add labels
plt.xlabel(r'sit( PIOMAS ) (m)')
plt.ylabel(r'sit( Submarine ) (m)')

### Adjust plot
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure + 'scattertest.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
###########################################################################

#fig = plt.figure()
#ax = plt.subplot(111)
#
#### Calculate trends
#varx = np.ravel(sitp2)
#vary = np.ravel(siti)
#
#mask = ~np.isnan(varx) & ~np.isnan(vary)
#slope, intercept, r_value, p_value, std_err = sts.linregress(varx[mask], vary[mask])
#line2 = slope*timex + intercept
#
#plt.plot(timex,timey,color='k',linewidth=2,zorder=1)
#plt.plot(line2,color='r',zorder=7)
#
#plt.scatter(sitp2[0,:,:],siti[0,:,:],label='2004',color='seagreen',zorder=2)
#plt.scatter(sitp2[1,:,:],siti[1,:,:],label='2005',color='rosybrown',zorder=3)
#plt.scatter(sitp2[2,:,:],siti[2,:,:],label='2006',color='sandybrown',zorder=4)
#plt.scatter(sitp2[3,:,:],siti[3,:,:],label='2007',color='steelblue',zorder=5)
#plt.scatter(sitp2[4,:,:],siti[4,:,:],label='2008',color='olive',zorder=6)
#plt.scatter(sitp2[5,:,:],siti[5,:,:],label='2009',color='Aquamarine',zorder=7)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#
#### Add legend
#plt.legend(shadow=False,fontsize=9,loc='center',
#                       fancybox=True,ncol=1,bbox_to_anchor=(0.93,0.14),
#                        frameon=False)
#                        
#plt.grid(color='darkgrey',linewidth=0.4)
#
#### Add limits
#plt.xlim([0,7])
#plt.ylim([0,7])
#plt.xticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
#plt.yticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
#
#### Add labels
##plt.xlabel(r'sit(PIOMAS) (m)')
##plt.ylabel(r'sit(ICESat-J) (m)')
#
#### Adjust plot
#fig.subplots_adjust(top=0.95)
#fig.subplots_adjust(bottom=0.15)
#
#plt.savefig(directoryfigure + 'scattertest2.png',dpi=300)

fig = plt.figure()
#ax = plt.subplot(111)

nullfmt = NullFormatter()

### definitions for the axes
left,width = 0.1,0.65
bottom,height = 0.1,0.65
bottom_h = left_h = left + width + 0.0

rect_scatter = [left,bottom,width,height]
rect_histx = [left,bottom_h,width,0.2]
rect_histy = [left_h,bottom,0.2,height]

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

### Calculate trends
varx = np.ravel(sitp2)
vary = np.ravel(siti)

#mask = ~np.isnan(varx) & ~np.isnan(vary)
mask = np.isfinite(varx) & np.isfinite(vary)
slope, intercept, r_value, p_value, std_err = sts.linregress(varx[mask], vary[mask])
line2 = slope*timex + intercept

fit = np.polyfit(varx[mask],vary[mask],2)
m = fit[0]
b = fit[1]
c = fit[2]

linetest = m*timex**2 + b*timex + c

axScatter.plot(timex,timey,color='k',linewidth=3,zorder=1)
#axScatter.plot(timex,linetest,color='r',zorder=7)

smoothed = sm.nonparametric.lowess(vary,varx)
axScatter.plot(smoothed[:,0],smoothed[:,1],color='r',zorder=8)

axScatter.scatter(sitp2[0,:,:],siti[0,:,:],label='2004',color='olive',zorder=2,s=5)
axScatter.scatter(sitp2[1,:,:],siti[1,:,:],label='2005',color='rosybrown',zorder=3,s=5)
axScatter.scatter(sitp2[2,:,:],siti[2,:,:],label='2006',color='sandybrown',zorder=4,s=5)
axScatter.scatter(sitp2[3,:,:],siti[3,:,:],label='2007',color='steelblue',zorder=5,s=5)
axScatter.scatter(sitp2[4,:,:],siti[4,:,:],label='2008',color='Aquamarine',zorder=6,s=5)
axScatter.scatter(sitp2[5,:,:],siti[5,:,:],label='2009',color='seagreen',zorder=7,s=5)

axScatter.set_xlim((0,7))
axScatter.set_ylim((0,7))
axScatter.set_xticklabels(map(str,np.arange(0,8,1)))
axScatter.set_yticklabels(map(str,np.arange(0,8,1)))

binwidth=0.25
bins = np.arange(0,7 + binwidth,binwidth)
n,bins,patches = axHistx.hist(varx[mask],bins=bins,normed=True,
             facecolor='dimgrey',edgecolor='w',alpha=1,linewidth=0.35)
n,binsy,patches = axHisty.hist(vary[mask],bins=bins,normed=True,
                               orientation='horizontal',
                               facecolor='dimgrey',edgecolor='w',
                               alpha=1,linewidth=0.35)

mu,sigma = sts.norm.fit(varx[mask])             
y = mlab.normpdf(bins,mu,sigma)
lx = axHistx.plot(bins,y,'k-',linewidth=0.5)

mu2,sigma2 = sts.norm.fit(vary[mask])
y2 = mlab.normpdf(binsy,mu2,sigma2)
ly = axHisty.plot(y2,binsy,'k-',linewidth=0.5)

axHistx.spines['top'].set_color('none')
axHistx.spines['right'].set_color('none')
axHistx.spines['left'].set_color('none')
axHistx.spines['bottom'].set_color('none')
axHisty.spines['top'].set_color('none')
axHisty.spines['right'].set_color('none')
axHisty.spines['left'].set_color('none')
axHisty.spines['bottom'].set_color('none')
axScatter.spines['top'].set_color('none')
axScatter.spines['right'].set_color('none')

axHistx.xaxis.set_ticks([])
axHistx.yaxis.set_ticks([])
axHisty.xaxis.set_ticks([])
axHisty.yaxis.set_ticks([])

axHistx.set_yticklabels([])
axHisty.set_xticklabels([])

axScatter.xaxis.set_ticks_position('bottom')
axScatter.yaxis.set_ticks_position('left')

#### Add legend
axScatter.legend(shadow=False,fontsize=7,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.93,0.18),
                        frameon=False)
                        
axScatter.grid(color='dimgrey',linewidth=0.4)

#### Add labels
axScatter.set_xlabel(r'\textbf{sit( PIOMAS )(m)}')
axScatter.set_ylabel(r'\textbf{sit( ICESat-J )(m)}')

#### Adjust plot
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.15)

### Add text
axHistx.text(5.05,0.05,r'*LOWESS Smoothing',fontsize=8)

plt.savefig(directoryfigure + 'scattertest2.png',dpi=800)

###########################################################################
###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
axb = fig.add_subplot(312)

datai = [varx[mask],vary[mask]]

vp = plt.violinplot(datai,showmeans=True,showmedians=False,vert=False,widths=0.6)

axb.spines['top'].set_color('none')
axb.spines['right'].set_color('none')
axb.spines['left'].set_color('none')
axb.spines['bottom'].set_color('none')
axb.xaxis.set_ticks_position('bottom')
axb.tick_params(left='off',right='off',bottom='off')
plt.setp(axb,xticks=[])

axb.set_aspect(1.9)

plt.setp(axb,yticks=[y+1 for y in range(len(datai))],
                     yticklabels=['PIOMAS','ICESat-J'])
                     
for i in vp['bodies']:
    i.set_edgecolor('darkgrey')  
vp['cbars'].set_color('k')
vp['cmaxes'].set_color('w')
vp['cmins'].set_color('w')
vp['cmeans'].set_color('k')
vp['cmeans'].set_linewidth(2)
vp['cmaxes'].set_linewidth(0.5)        
vp['cmins'].set_linewidth(0.5)       
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')          
vp['bodies'][0].set_facecolor('seagreen')
vp['bodies'][1].set_facecolor('goldenrod')                      

plt.xlim([0,9]) 
plt.xticks(np.arange(0,10,1),[]) 

###########################################################################
###########################################################################

axb = fig.add_subplot(311)

varxb = np.ravel(sitp)
varyb = np.ravel(sitb)
maskb = np.isfinite(varxb) & np.isfinite(varyb)

datab=[varxb[maskb],varyb[maskb]]

vp = plt.violinplot(datab,showmeans=True,showmedians=False,vert=False,widths=0.6)

axb.spines['top'].set_color('none')
axb.spines['right'].set_color('none')
axb.spines['left'].set_color('none')
axb.spines['bottom'].set_color('none')
axb.xaxis.set_ticks_position('bottom')
axb.tick_params(left='off',right='off',bottom='off')

axb.set_aspect(1.9)

plt.setp(axb,yticks=[y+1 for y in range(len(datai))],
                     yticklabels=['PIOMAS','Submarine'])
                     
for i in vp['bodies']:
    i.set_edgecolor('darkgrey')  
vp['cbars'].set_color('k')
vp['cmaxes'].set_color('w')
vp['cmins'].set_color('w')
vp['cmeans'].set_color('k')
vp['cmaxes'].set_linewidth(0.5)        
vp['cmins'].set_linewidth(0.5) 
vp['cmeans'].set_linewidth(2)
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')              
vp['bodies'][0].set_facecolor('seagreen')
vp['bodies'][1].set_facecolor('darkblue')                      
 
plt.xlim([0,9])                    
plt.xticks(np.arange(0,10,1),[])      
axb.tick_params('both',length=5.5,width=1,which='major')  

###########################################################################
###########################################################################

axb = fig.add_subplot(313)

varxc = np.ravel(sitp3)
varyc = np.ravel(sitc)
maskc = np.isfinite(varxc) & np.isfinite(varyc)

datac=[varxc[maskc],varyc[maskc]]

vp = plt.violinplot(datac,showmeans=True,showmedians=False,vert=False,widths=0.6)

axb.spines['top'].set_color('none')
axb.spines['right'].set_color('none')
axb.spines['left'].set_color('none')
axb.xaxis.set_ticks_position('bottom')
axb.tick_params(left='off',right='off')

axb.set_aspect(1.9)

plt.setp(axb,yticks=[y+1 for y in range(len(datai))],
                     yticklabels=['PIOMAS','CryoSat-2'])
                     
for i in vp['bodies']:
    i.set_edgecolor('darkgrey')  
vp['cbars'].set_color('k')
vp['cmaxes'].set_color('w')
vp['cmins'].set_color('w')
vp['cmeans'].set_color('k')
vp['cmaxes'].set_linewidth(0.5)        
vp['cmins'].set_linewidth(0.5) 
vp['cmeans'].set_linewidth(2)
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')              
vp['bodies'][0].set_facecolor('seagreen')
vp['bodies'][1].set_facecolor('darkred')                      
                                         
plt.xticks(np.arange(0,10,1),map(str,np.arange(0,10,1)))    
plt.xlabel(r'Thickness (m)',fontsize=10)   
axb.tick_params('both',length=5.5,width=1,which='major')

###########################################################################
###########################################################################

masking = meansitb.copy()
masking[np.where(masking == 0.)] = np.nan
masking[np.where(np.isfinite(masking))] = 1.

a2 = plt.axes([.67, .68, .29, .22], axisbg='w')   
m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'dimgrey',linewidth=0.2)
m.drawlsmask(land_color='dimgrey',ocean_color='snow')
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.15)
m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.15,
                fontsize=3)

cs = m.contourf(lon,lat,masking,np.arange(0,3,1),
                latlon=True,colors='darkblue',alpha=0.5)
                
masking = np.nanmean(siti,axis=0)
masking[np.where(masking == 0.)] = np.nan
masking[np.where(np.isfinite(masking))] = 1.

a2 = plt.axes([.67, .41, .29, .22], axisbg='w')   
m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'dimgrey',linewidth=0.2)
m.drawlsmask(land_color='dimgrey',ocean_color='snow')
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.15)
m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.15,
                fontsize=3)

cs = m.contourf(lon,lat,masking,np.arange(0,3,1),
                latlon=True,colors='goldenrod',alpha=0.5)     
                
masking = np.nanmean(sitc,axis=0)
masking[np.where(masking == 0.)] = np.nan
masking[np.where(np.isfinite(masking))] = 1.

a2 = plt.axes([.67, .14, .29, .22], axisbg='w')   
m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'dimgrey',linewidth=0.2)
m.drawlsmask(land_color='dimgrey',ocean_color='snow')
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.15)
m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.15,
                fontsize=3)

cs = m.contourf(lon,lat,masking,np.arange(0,3,1),
                latlon=True,colors='darkred',alpha=0.5)                 

plt.savefig(directoryfigure + 'boxtest2.png',dpi=800)

###########################################################################
###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
#ax = plt.subplot(111)

nullfmt = NullFormatter()

### definitions for the axes
left,width = 0.1,0.65
bottom,height = 0.1,0.65
bottom_h = left_h = left + width + 0.0

rect_scatter = [left,bottom,width,height]
rect_histx = [left,bottom_h,width,0.2]
rect_histy = [left_h,bottom,0.2,height]

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

### Calculate trends
varx = np.ravel(sitp3)
vary = np.ravel(sitc)

#mask = ~np.isnan(varx) & ~np.isnan(vary)
mask = np.isfinite(varx) & np.isfinite(vary)
slope, intercept, r_value, p_value, std_err = sts.linregress(varx[mask], vary[mask])
line2 = slope*timex + intercept

fit = np.polyfit(varx[mask],vary[mask],1)
m = fit[0]
b = fit[1]

linetest = m*timex+b

smoothed = sm.nonparametric.lowess(vary,varx)
axScatter.plot(smoothed[:,0],smoothed[:,1],color='r',zorder=8)

axScatter.plot(timex,timey,color='k',linewidth=3,zorder=1)
#axScatter.plot(timex,linetest,color='r',zorder=7)

axScatter.scatter(sitp3[0,:,:],sitc[0,:,:],label='2011',color='olive',zorder=2,s=5)
axScatter.scatter(sitp3[1,:,:],sitc[1,:,:],label='2012',color='rosybrown',zorder=3,s=5)
axScatter.scatter(sitp3[2,:,:],sitc[2,:,:],label='2013',color='sandybrown',zorder=4,s=5)
axScatter.scatter(sitp3[3,:,:],sitc[3,:,:],label='2014',color='steelblue',zorder=5,s=5)
axScatter.scatter(sitp3[4,:,:],sitc[4,:,:],label='2015',color='Aquamarine',zorder=6,s=5)

axScatter.set_xlim((0,7))
axScatter.set_ylim((0,7))
axScatter.set_xticklabels(map(str,np.arange(0,8,1)))
axScatter.set_yticklabels(map(str,np.arange(0,8,1)))

binwidth=0.25
bins = np.arange(0,7 + binwidth,binwidth)
n,bins,patches = axHistx.hist(varx[mask],bins=bins,normed=True,
             facecolor='dimgrey',edgecolor='w',alpha=1,linewidth=0.35)
n,binsy,patches = axHisty.hist(vary[mask],bins=bins,normed=True,
                               orientation='horizontal',
                               facecolor='dimgrey',edgecolor='w',
                               alpha=1,linewidth=0.35)

mu,sigma = sts.norm.fit(varx[mask])             
y = mlab.normpdf(bins,mu,sigma)
lx = axHistx.plot(bins,y,'k-',linewidth=0.5)

mu2,sigma2 = sts.norm.fit(vary[mask])
y2 = mlab.normpdf(binsy,mu2,sigma2)
ly = axHisty.plot(y2,binsy,'k-',linewidth=0.5)

axHistx.spines['top'].set_color('none')
axHistx.spines['right'].set_color('none')
axHistx.spines['left'].set_color('none')
axHistx.spines['bottom'].set_color('none')
axHisty.spines['top'].set_color('none')
axHisty.spines['right'].set_color('none')
axHisty.spines['left'].set_color('none')
axHisty.spines['bottom'].set_color('none')
axScatter.spines['top'].set_color('none')
axScatter.spines['right'].set_color('none')

axHistx.xaxis.set_ticks([])
axHistx.yaxis.set_ticks([])
axHisty.xaxis.set_ticks([])
axHisty.yaxis.set_ticks([])

axHistx.set_yticklabels([])
axHisty.set_xticklabels([])

axScatter.xaxis.set_ticks_position('bottom')
axScatter.yaxis.set_ticks_position('left')

#### Add legend
axScatter.legend(shadow=False,fontsize=7,loc='center',
                       fancybox=True,ncol=1,bbox_to_anchor=(0.93,0.18),
                        frameon=False)
                        
axScatter.grid(color='dimgrey',linewidth=0.4)

### Add text
axHistx.text(5.05,0.05,r'*LOWESS Smoothing',fontsize=8)

#### Add labels
axScatter.set_xlabel(r'\textbf{sit( PIOMAS )(m)}')
axScatter.set_ylabel(r'\textbf{sit( CryoSat-2 )(m)}')
#
#### Adjust plot
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure + 'scattertest3.png',dpi=800)

print 'Completed: Script done!'