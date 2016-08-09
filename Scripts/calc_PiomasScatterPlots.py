"""
Script uses PIOMAS data for snow cover, oflux, sea ice thickness,
and sea ice concentration to make scatter plots
 
Author : Zachary Labe
Date : 15 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as c
import datetime

### Define directories
directorydatasit = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directorydatasic = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/SeaIceConcentration/' 
directorydatasnow = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/SnowCover/'  
directorydataoflux = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/OceanFlux/' 
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

print '\n' '--- PIOMAS Scatter Plots (%s) ---' '\n' % titletime 

def piomasReader(directory,var,year,month,time):
    """
    Reads piomas data for sea ice thickness over 1979-2015
    """
    
    ### Enter filename
    if var == 'thick':
        filename = 'piomas_regrid_sit_19792015.nc'   
    elif var == 'sic':
        filename = 'piomas_regrid_sic_19792015.nc'
    elif var == 'snow':
        filename = 'piomas_regrid_snow_19792015.nc'
    elif var == 'oflux':
        filename = 'piomas_regrid_oflux_19792004.nc'

    ### Retrieve data
    data = Dataset(directory + filename)
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    
    if time == True:
        if var == 'thick':
            vari = data.variables['newthickness'][year,month,:,:]   
        elif var == 'sic':
            vari = data.variables['sic'][year,month,:,:]
        elif var == 'snow':
            vari = data.variables['snow'][year,month,:,:]            
        elif var == 'oflux':
            vari = data.variables['oflux'][year,month,:,:]
    else:
        if var == 'thick':
            vari = data.variables['newthickness'][:,:,:,:]   
        elif var == 'sic':
            vari = data.variables['sic'][:,:,:,:]
        elif var == 'snow':
            vari = data.variables['snow'][:,:,:,:]            
        elif var == 'oflux':
            vari = data.variables['oflux'][:,:,:,:]
        
            ### Convert to new units (seconds --> days)
            vari = vari * 60 * 60 * 24
        
    data.close()
    
    print 'Completed: PIOMAS data read!'
    return lat,lon,vari
    
def scatterPlot(sit,sic,snow,flux,varx,vary,year,month):

    fig = plt.figure()
    ax = plt.subplot(111)
    
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ### Define colormap
    def colormapSIT():
        cmap1 = plt.get_cmap('BuPu')
        cmap2 = plt.get_cmap('RdPu_r')
        cmap3 = plt.get_cmap('gist_heat_r')
        cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
        cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
        cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
        cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
        return cms_sit
        
    def colormapSIC():
        cmap = plt.get_cmap('RdPu')
        cmaplist = [cmap(i) for i in xrange(0,cmap.N-20)]
        cms_sic = c.ListedColormap(cmaplist)
        return cms_sic
        
    def colormapOflux():
        cmap = plt.get_cmap('hot_r')
        cmaplist = [cmap(i) for i in xrange(60,cmap.N-40)]
        cms_sic = c.ListedColormap(cmaplist)
        return cms_sic
        
    def colormapSnow():
        cmap1 = plt.get_cmap('viridis')
        cms_snow = cmap1 
        return cms_snow
    
    if varx == 'thick':
        x = sit
        xlabel = r'\textbf{Sea Ice Thickness (meters)}'
        cm = colormapSIT()
        lim = np.arange(0,8,1)
        limx1 = 0
        limx2 = 7
    elif varx == 'sic':
        x = sic
        xlabel = r'\textbf{Sea Ice Concentration (fraction \%)}'
        cm = colormapSIC()
        lim = np.arange(0,1.1,0.25)
        limx1 = 0
        limx2 = 1
    elif varx == 'snow':
        x = snow
        xlabel = r'\textbf{Snow Cover (meters)}'
        cm = colormapSnow()
        lim = np.arange(0,0.5,0.1)
        limx1 = 0
        limx2 = 0.4
    elif varx == 'oflux':
        x = flux 
        xlabel = r'\textbf{Ocean Heat Flux (meters of ice/day)}'
        cm = colormapOflux()
        lim = np.arange(0,0.2,0.2)
        limx1 = 0
        limx2 = 0.14
        
    if vary == 'thick':
        y = sit
        ylabel = r'\textbf{Sea Ice Thickness (meters)}'
        limy1 = 0
        limy2 = 7
    elif vary == 'sic':
        y = sic
        ylabel = r'\textbf{Sea Ice Concentration (fraction \%)}'
        limy1 = 0
        limy2 = 1
    elif vary == 'snow':
        y = snow
        ylabel = r'\textbf{Snow Cover (meters)}'
        limy1 = 0
        limy2 = 0.4
    elif vary == 'oflux':
        y = flux
        ylabel = r'\textbf{Ocean Heat Flux (meters of ice/day)}'
        limy1 = 0
        limy2 = 0.14
        
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
    
    ### Adjust selected year and month
    x = x[year,month,:,:]
    y = y[year,month,:,:]
    
    ### Log/ln adjustments
    y = np.log10(y)
#    y = np.log(y)
    
    cs = plt.scatter(x[:,:],y[:,:],c=x,cmap=cm)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([limx1,limx2])
#    plt.ylim([limy1,limy2])
    
    cbar = plt.colorbar(cs,pad=.01)
    cbar.set_ticks(lim)
    cbar.set_ticklabels(lim)
    cbar.set_label(xlabel)
    
    plt.title(r'\textbf{PIOMAS %s %s}' % (datemo,dateyr),fontsize=18)
    
#    plt.savefig(directoryfigure + '/ScatterPlots/piomas_%s%s_%s%s_log.png' % (varx,vary,datemo[0:3],str(dateyr)[-2:]),
#                dpi=300)
    
    print 'Completed: Scatter plot!'
    
def meanScatter(sic,sit,snow,flux,year,years,month):
    meansit = np.nanmean(np.nanmean(sit,axis=2),axis=2)
    meansic = np.nanmean(np.nanmean(sic,axis=2),axis=2)
    meansityr = np.nanmean(meansit,axis=1)
    meansicyr = np.nanmean(meansic,axis=1)
    
    meansnow = np.nanmean(np.nanmean(snow,axis=2),axis=2)
    meanflux = np.nanmean(np.nanmean(flux,axis=2),axis=2)
    meansnowyr = np.nanmean(meansnow,axis=1)
    meanfluxyr = np.nanmean(meanflux,axis=1)
    
    fig = plt.figure()
    
    dateyr = years[year]    
    datemo = datetime.date(dateyr,month+1,1).strftime('%B')
    
    ### Call parameters
    plt.rcParams['text.usetex']=True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = 'Avant Garde'
    
    ax = plt.subplot(221)
    plt.scatter(meansityr,meansicyr,color='darkviolet',edgecolor='k')
    plt.xlim([0.6,1.6])
    plt.ylabel(r'\textbf{SIC}')
    
    ax = plt.subplot(222)
    plt.scatter(meansityr[:26],meansnowyr,color='darkslategrey',edgecolor='k')
    plt.xlim([0.6,1.6])
    plt.ylabel(r'\textbf{Snow Cover}')
    plt.xlabel(r'\textbf{SIT}')
    
    ax = plt.subplot(223)
    plt.scatter(meansityr[:26],meanfluxyr,color='darkorange',edgecolor='k')
    plt.xlim([0.6,1.6])
    plt.ylabel(r'\textbf{Ocean Flux}')
    plt.xlabel(r'\textbf{SIT}')
    
    fig.subplots_adjust(wspace=0.3)
    fig.subplots_adjust(hspace=0.3)
    fig.suptitle(r'\textbf{PIOMAS Temporal Relationship with SIT}',
                 fontsize=15)
    
    plt.savefig(directoryfigure + '/ScatterPlots/subplot_scatter.png',
                dpi = 300)
                
    print 'Completed: Temporal subplot!'
    
### Call functions
year = 0   
month = 2
varx = 'thick'
vary = 'oflux'
time = False
    
lat,lon,sit = piomasReader(directorydatasit,'thick',year,month,time)
lat,lon,sic = piomasReader(directorydatasic,'sic',year,month,time)
lat,lon,snow = piomasReader(directorydatasnow,'snow',year,month,time)
lat,lon,flux = piomasReader(directorydataoflux,'oflux',year,month,time)
scatterPlot(sit[:26],sic,snow,flux,varx,vary,year,month)
#meanScatter(sic,sit,snow,flux,year,years,month)

print 'Completed: Script finished!'
