"""
Script creates scatter plots between March data for SIT using PIOMAS
and ICESat-J and CryoSat-2. Subplots are for each individual year 
over the 2004-2015 time frame.
 
Author : Zachary Labe
Date : 8 August 2016
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

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'  
directorydata2= '/home/zlabe/Surtsey/seaice_obs/PIOMAS/SeaIceConcentration/' 
#directoryfigure = '/home/zlabe/Desktop/'
directoryfigure = '/home/zlabe/Documents/Research/SeaIceThickness/Figures/Testsq/' 

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

print '\n' '--- PIOMAS/Satellites Yearly Scatter Plots (%s) ---' '\n' % titletime 

def piomasReader(directory,directory2,segment,years):
    
    filename = 'piomas_regrid_March_19792015.nc'
    
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    
    filename2 = 'piomas_regrid_sic_19792015.nc'
    data2 = Dataset(directory2 + filename2)
    
    if segment == 'sub':    # 1986-1994  
        timeslice = np.where((years >= 1986) & (years <= 1994))[0]
        sitp = data2.variables['thick'][timeslice,:,:] 
        sicp = data2.variables['sic'][timeslice,2,:,:] 
        sitp = sitp / sicp
    elif segment == 'icej':    # 2004-2009
        timeslice = np.where((years >= 2004) & (years <= 2009))[0]
        sitp = data.variables['thick'][timeslice,:,:]
        sicp = data2.variables['sic'][timeslice,2,:,:]
        sitp = sitp / sicp
    elif segment == 'cryo':    # 2011-2015
        timeslice = np.where((years >= 2011) & (years <= 2015))[0]
        sitp = data.variables['thick'][timeslice,:,:]
        sicp = data2.variables['sic'][timeslice,2,:,:]
        sitp = sitp / sicp
    else:
        sitp = data.variables['thick'][:,:,:]
        sicp = data2.variables['sic'][:,2,:,:]
        sitp = sitp / sicp
        
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
lat,lon,sitpi = piomasReader(directorydata,directorydata2,'icej',years)
lat,lon,sitpc = piomasReader(directorydata,directorydata2,'cryo',years)
lat,lon,siti,sitc = icesatReader(directorydata)

timex = np.arange(0,8,1)
timey = np.arange(0,8,1)

### Calculate trends
varx1 = np.ravel(sitpi[0,:,:])
vary1 = np.ravel(siti[0,:,:])
mask1 = np.isfinite(varx1) & np.isfinite(vary1)
fit1 = np.polyfit(varx1[mask1],vary1[mask1],1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx1[mask1],
                                                          vary1[mask1])                                                        
m1 = fit1[0]
b1 = fit1[1]
line1 = m1*timex+b1

varx2 = np.ravel(sitpi[1,:,:])
vary2 = np.ravel(siti[1,:,:])
mask2 = np.isfinite(varx2) & np.isfinite(vary2)
fit2 = np.polyfit(varx2[mask2],vary2[mask2],1)
slope,intercept,r2,p_value,std_err = sts.stats.linregress(varx2[mask2],
                                                          vary2[mask2])
m2 = fit2[0]
b2 = fit2[1]
line2 = m2*timex+b2

varx3 = np.ravel(sitpi[2,:,:])
vary3 = np.ravel(siti[2,:,:])
mask3 = np.isfinite(varx3) & np.isfinite(vary3)
fit3 = np.polyfit(varx3[mask3],vary3[mask3],1)
slope,intercept,r3,p_value,std_err = sts.stats.linregress(varx3[mask3],
                                                          vary3[mask3])
m3 = fit3[0]
b3 = fit3[1]
line3 = m3*timex+b3

varx4 = np.ravel(sitpi[3,:,:])
vary4 = np.ravel(siti[3,:,:])
mask4 = np.isfinite(varx4) & np.isfinite(vary4)
fit4 = np.polyfit(varx4[mask4],vary4[mask4],1)
slope,intercept,r4,p_value,std_err = sts.stats.linregress(varx4[mask4],
                                                          vary4[mask4])
m4 = fit4[0]
b4 = fit4[1]
line4 = m4*timex+b4

varx5 = np.ravel(sitpi[4,:,:])
vary5 = np.ravel(siti[4,:,:])
mask5 = np.isfinite(varx5) & np.isfinite(vary5)
fit5 = np.polyfit(varx5[mask5],vary5[mask5],1)
slope,intercept,r5,p_value,std_err = sts.stats.linregress(varx5[mask5],
                                                          vary5[mask5])
m5 = fit5[0]
b5 = fit5[1]
line5 = m5*timex+b5

varx6 = np.ravel(sitpi[5,:,:])
vary6 = np.ravel(siti[5,:,:])
mask6 = np.isfinite(varx6) & np.isfinite(vary6)
fit6 = np.polyfit(varx6[mask6],vary6[mask6],1)
slope,intercept,r6,p_value,std_err = sts.stats.linregress(varx6[mask6],
                                                          vary6[mask6])
m6 = fit6[0]
b6 = fit6[1]
line6 = m6*timex+b6

### ICESat-J yearly plots
fig = plt.figure()

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

ax1 = plt.subplot(231)
plt.scatter(sitpi[0,:,:],siti[0,:,:],label='2004',color='steelblue',zorder=2,s=2)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax1.plot(timex,line1,color='r',zorder=7)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'2004',color='steelblue',fontsize=11)
ax1.text(4.7,0.5,r'r$^2$= %s' % round(r1**2,2),color='k',fontsize=7)

ax2 = plt.subplot(232)
plt.scatter(sitpi[1,:,:],siti[1,:,:],label='2005',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.plot(timex,line2,color='r',zorder=7)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'2005',color='steelblue',fontsize=11)
ax2.text(4.7,0.5,r'r$^2$= %s' % round(r2**2,2),color='k',fontsize=7)

ax3 = plt.subplot(233)
plt.scatter(sitpi[2,:,:],siti[2,:,:],label='2006',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,line3,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'2006',color='steelblue',fontsize=11)
ax3.text(4.7,0.5,r'r$^2$= %s' % round(r3**2,2),color='k',fontsize=7)

ax4 = plt.subplot(234)
plt.scatter(sitpi[3,:,:],siti[3,:,:],label='2007',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax4.plot(timex,line4,color='r',zorder=7)
ax4.set_xlim((0,7))
ax4.set_ylim((0,7))
ax4.set_xticklabels(map(str,np.arange(0,8,1)))
ax4.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax4, ['left', 'bottom'])
ax4.spines['top'].set_color('none')
ax4.spines['right'].set_color('none')
ax4.text(0,6.2,r'2007',color='steelblue',fontsize=11)
ax4.text(4.7,0.5,r'r$^2$= %s' % round(r4**2,2),color='k',fontsize=7)

ax5 = plt.subplot(235)
plt.scatter(sitpi[4,:,:],siti[4,:,:],label='2008',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax5.plot(timex,line5,color='r',zorder=7)
ax5.set_xlim((0,7))
ax5.set_ylim((0,7))
ax5.set_xticklabels(map(str,np.arange(0,8,1)))
ax5.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax5, ['left', 'bottom'])
ax5.spines['top'].set_color('none')
ax5.spines['right'].set_color('none')
ax5.text(0,6.2,r'2008',color='steelblue',fontsize=11)
ax5.text(4.7,0.5,r'r$^2$= %s' % round(r5**2,2),color='k',fontsize=7)

ax6 = plt.subplot(236)
plt.scatter(sitpi[5,:,:],siti[5,:,:],label='2009',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax6.plot(timex,line6,color='r',zorder=7)
ax6.set_xlim((0,7))
ax6.set_ylim((0,7))
ax6.set_xticklabels(map(str,np.arange(0,8,1)))
ax6.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax6, ['left', 'bottom'])
ax6.spines['top'].set_color('none')
ax6.spines['right'].set_color('none')
ax6.text(0,6.2,r'2009',color='steelblue',fontsize=11)
ax6.text(4.7,0.5,r'r$^2$= %s' % round(r6**2,2),color='k',fontsize=7)

fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-9.4,-3.2,r'\textbf{sit( PIOMAS )(m)}',fontsize=16,rotation=0)
plt.text(-21,11.5,r'\textbf{sit( ICESat-J )(m)}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'yrtest.png',dpi=600)

###########################################################################
###########################################################################
###########################################################################

### Calculate trends
varx1 = np.ravel(sitpc[0,:,:])
vary1 = np.ravel(sitc[0,:,:])
mask1 = np.isfinite(varx1) & np.isfinite(vary1)
fit1 = np.polyfit(varx1[mask1],vary1[mask1],1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx1[mask1],
                                                          vary1[mask1]) 
m1 = fit1[0]
b1 = fit1[1]
line1 = m1*timex+b1

varx2 = np.ravel(sitpc[1,:,:])
vary2 = np.ravel(sitc[1,:,:])
mask2 = np.isfinite(varx2) & np.isfinite(vary2)
fit2 = np.polyfit(varx2[mask2],vary2[mask2],1)
slope,intercept,r2,p_value,std_err = sts.stats.linregress(varx2[mask2],
                                                          vary2[mask2]) 
m2 = fit2[0]
b2 = fit2[1]
line2 = m2*timex+b2

varx3 = np.ravel(sitpc[2,:,:])
vary3 = np.ravel(sitc[2,:,:])
mask3 = np.isfinite(varx3) & np.isfinite(vary3)
fit3 = np.polyfit(varx3[mask3],vary3[mask3],1)
slope,intercept,r3,p_value,std_err = sts.stats.linregress(varx3[mask3],
                                                          vary3[mask3]) 
m3 = fit3[0]
b3 = fit3[1]
line3 = m3*timex+b3

varx4 = np.ravel(sitpc[3,:,:])
vary4 = np.ravel(sitc[3,:,:])
mask4 = np.isfinite(varx4) & np.isfinite(vary4)
fit4 = np.polyfit(varx4[mask4],vary4[mask4],1)
slope,intercept,r4,p_value,std_err = sts.stats.linregress(varx4[mask4],
                                                          vary4[mask4]) 
m4 = fit4[0]
b4 = fit4[1]
line4 = m4*timex+b4

varx5 = np.ravel(sitpc[4,:,:])
vary5 = np.ravel(sitc[4,:,:])
mask5 = np.isfinite(varx5) & np.isfinite(vary5)
fit5 = np.polyfit(varx5[mask5],vary5[mask5],1)
slope,intercept,r5,p_value,std_err = sts.stats.linregress(varx5[mask5],
                                                          vary5[mask5]) 
m5 = fit5[0]
b5 = fit5[1]
line5 = m5*timex+b5

### CryoSat-2 yearly plots
fig = plt.figure()

ax1 = plt.subplot(231)
plt.scatter(sitpc[0,:,:],sitc[0,:,:],label='2011',color='steelblue',zorder=2,s=2)
ax1.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax1.plot(timex,line1,color='r',zorder=7)
ax1.set_xlim((0,7))
ax1.set_ylim((0,7))
ax1.set_xticklabels(map(str,np.arange(0,8,1)))
ax1.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax1, ['left', 'bottom'])
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.text(0,6.2,r'2011',color='steelblue',fontsize=11)
ax1.text(4.7,0.5,r'r$^2$= %s' % round(r1**2,2),color='k',fontsize=7)

ax2 = plt.subplot(232)
plt.scatter(sitpc[1,:,:],sitc[1,:,:],label='2012',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax2.plot(timex,line2,color='r',zorder=7)
ax2.set_xlim((0,7))
ax2.set_ylim((0,7))
ax2.set_xticklabels(map(str,np.arange(0,8,1)))
ax2.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax2, ['left', 'bottom'])
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.text(0,6.2,r'2012',color='steelblue',fontsize=11)
ax2.text(4.7,0.5,r'r$^2$= %s' % round(r2**2,2),color='k',fontsize=7)

ax3 = plt.subplot(233)
plt.scatter(sitpc[2,:,:],sitc[2,:,:],label='2013',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax3.plot(timex,line3,color='r',zorder=7)
ax3.set_xlim((0,7))
ax3.set_ylim((0,7))
ax3.set_xticklabels(map(str,np.arange(0,8,1)))
ax3.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax3, ['left', 'bottom'])
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.text(0,6.2,r'2013',color='steelblue',fontsize=11)
ax3.text(4.7,0.5,r'r$^2$= %s' % round(r3**2,2),color='k',fontsize=7)

ax4 = plt.subplot(234)
plt.scatter(sitpc[3,:,:],sitc[3,:,:],label='2014',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax4.plot(timex,line4,color='r',zorder=7)
ax4.set_xlim((0,7))
ax4.set_ylim((0,7))
ax4.set_xticklabels(map(str,np.arange(0,8,1)))
ax4.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax4, ['left', 'bottom'])
ax4.spines['top'].set_color('none')
ax4.spines['right'].set_color('none')
ax4.text(0,6.2,r'2014',color='steelblue',fontsize=11)
ax4.text(4.7,0.5,r'r$^2$= %s' % round(r4**2,2),color='k',fontsize=7)

ax5 = plt.subplot(235)
plt.scatter(sitpc[4,:,:],sitc[4,:,:],label='2015',color='steelblue',zorder=2,s=2)
plt.plot(timex,timey,color='k',linewidth=3,zorder=1)
ax5.plot(timex,line5,color='r',zorder=7)
ax5.set_xlim((0,7))
ax5.set_ylim((0,7))
ax5.set_xticklabels(map(str,np.arange(0,8,1)))
ax5.set_yticklabels(map(str,np.arange(0,8,1)))
adjust_spines(ax5, ['left', 'bottom'])
ax5.spines['top'].set_color('none')
ax5.spines['right'].set_color('none')
ax5.text(0,6.2,r'2015',color='steelblue',fontsize=11)
ax5.text(4.7,0.5,r'r$^2$= %s' % round(r5**2,2),color='k',fontsize=7)

fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(hspace=.4)
fig.subplots_adjust(bottom=0.15)

plt.text(-.2,-3.2,r'\textbf{sit( PIOMAS )}',fontsize=16,rotation=0)
plt.text(-12,11.5,r'\textbf{sit( CryoSat-2 )}',fontsize=16,rotation=90)

plt.savefig(directoryfigure + 'yrtest2.png',dpi=300)