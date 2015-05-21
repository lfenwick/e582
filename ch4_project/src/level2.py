from __future__ import division,print_function
import glob
import h5py
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib import pyplot as plt
import site
site.addsitedir('../utilities')
import reproject
from stopwatch import timeit
import compat
if compat.is_py3:
    import importlib
    reload=importlib.reload
reload(reproject)
from reproject import reproj_L1B,reproj_numba
import numpy as np
import numpy.ma as ma



def make_plot(lcc_values):
    """
      set up the basic map projection details with coastlines and meridians
      return the projection object for further plotting
    """
    proj = Basemap(**lcc_values)
    parallels = np.arange(-90, 90, 5)
    meridians = np.arange(0, 360, 5)
    proj.drawparallels(parallels, labels=[1, 0, 0, 0],
                       fontsize=10, latmax=90)
    proj.drawmeridians(meridians, labels=[0, 0, 0, 1],
                       fontsize=10, latmax=90)
    # draw coast & fill continents
    # map.fillcontinents(color=[0.25, 0.25, 0.25], lake_color=None) # coral
    proj.drawcoastlines(linewidth=1.5, linestyle='solid', color='k')
    return proj


def find_corners(lons, lats):
    """
      guess values for the upper right and lower left corners of the
      lat/lon grid and the grid center based on max/min lat lon in the
      data and return a dictionary that can be passed to Basemap to set
      the lcc projection.  Also return the smallest lat and lon differences
      to get a feeling for the image resolution
    """
    min_lat, min_lon = np.min(lats), np.min(lons)
    max_lat, max_lon = np.max(lats), np.max(lons)
    llcrnrlon, llcrnrlat = min_lon, min_lat
    urcrnrlon, urcrnrlat = max_lon, max_lat
    lon_res=np.min(np.abs(np.diff(lons.flat)))
    lat_res=np.min(np.abs(np.diff(lats.flat)))
    out=dict(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
             lat_1=llcrnrlat,lat_2=urcrnrlat,lat_0=(llcrnrlat+urcrnrlat)/2.,
             lon_0=(llcrnrlon + urcrnrlon)/2.)
    return(out,lon_res,lat_res)
    

if __name__ == "__main__":
    filename='A2010130213500.h5'
    path=glob.glob('../datasets/{}'.format(filename))[0]


    with h5py.File(path,'r') as f:
        lats=f['Navigation Data/latitude'][...]
        lons=f['Navigation Data/longitude'][...]
        chlor=f['Geophysical Data/chlor_a']
        chlor_array=chlor[...]
        chlor_bad_value=chlor.attrs['bad_value_scaled']

    #
    # take the last 500 rows and columns for testing
    #
    lim= -500
    the_slice=slice(lim,None)
    small_lons=lons[the_slice]
    small_lats=lats[the_slice]
    chlor_array=chlor_array[the_slice]
    
    plt.close('all')
    #
    # pixels without map projection or binning
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    out=ax.plot(small_lons,small_lats,'b+')
    ax.set_title('raw pixel center lat/lons')

    #
    # get the basic Basemap option dictionary
    #
    lcc_values,lon_res,lat_res=find_corners(small_lons,small_lats)
    lcc_values['resolution']='l'
    lcc_values['projection']='lcc'
    fig.savefig('nobins.png')
    
    #
    #pixels with map projection
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    #
    # tell Basemap what axis to plot into
    #
    lcc_values['ax']=ax
    proj=make_plot(lcc_values)
    x,y=proj(small_lons,small_lats)
    proj.ax.plot(x,y,'b+')
    proj.ax.set_title('pixel centers on lcc projection')
    proj.ax.figure.canvas.draw()
    fig.savefig('centers.png')
    
    cmap=cm.YlGn  #see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
    cmap.set_over('r')
    cmap.set_under('b')
    cmap.set_bad('0.75') #75% grey
    vmin= -1.5
    vmax= 1.5
    res=0.05
    the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
    #
    # gridded chlorophyll
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    lcc_values['ax']=ax
    proj=make_plot(lcc_values)

    lonlim=[np.min(lons), np.max(lons)]
    latlim=[np.min(lats), np.max(lats)]
    with timeit('chlor_slow'):
        chlor_a, longitude, latitude, bin_count = reproj_L1B(chlor_array,chlor_bad_value, small_lons, small_lats, lonlim, latlim, res)
    with timeit('chlor_numba'):
        chlor_a, longitude, latitude, bin_count = reproj_numba(chlor_array,chlor_bad_value, small_lons, small_lats, lonlim, latlim, res)
    with timeit('chlor_numba 2'):
        chlor_a, longitude, latitude, bin_count = reproj_numba(chlor_array,chlor_bad_value, small_lons, small_lats, lonlim, latlim, res)
    log_chlor=np.log10(chlor_a)
    mask=np.isnan(log_chlor)
    log_chlor=ma.array(log_chlor,mask=mask)
    x,y=proj(longitude,latitude)
    CS=proj.pcolormesh(x, y, log_chlor, cmap=cmap,norm=the_norm)
    CBar=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
    CBar.set_label('log10 (mg/m^3) of the chlorophyll a concentration')
    proj.ax.set_title('binned chlorophyll values on a regular grid')
    proj.ax.figure.canvas.draw()
    fig.savefig('chloro.png')

        
    #plt.show()





