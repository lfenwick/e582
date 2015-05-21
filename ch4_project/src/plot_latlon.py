from __future__ import division, print_function
from mpl_toolkits.basemap import Basemap
import h5py
from matplotlib import pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import cm
from matplotlib.colors import Normalize
import seaborn as sns

if __name__ == '__main__':

    in_file='navigated.h5'
    with h5py.File(in_file, "r") as f:
        lat = f['lattitude'][...]
        lon = f['longitude'][...]
        chlor=f['chlor_a_mean'][...]
        units=f['chlor_a_mean'].attrs['units']
        title=f['chlor_a_mean'].attrs['title']
        start=f.attrs['start_date']
        stop=f.attrs['end_date']

    lcc_values=dict(resolution='l',projection='lcc',
                    lat_1=30,lat_2=50,lat_0=45,lon_0=-135,
                    llcrnrlon=-133,llcrnrlat=40,
                    urcrnrlon=-118,urcrnrlat=50)
    proj=Basemap(**lcc_values)
    #
    # set up a yellow-green colormap
    # with reserved colors for over (red), under (blue) and missing (grey)
    # values
    #
    cmap=cm.YlGn  #see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
    cmap.set_over('r')
    cmap.set_under('b')
    cmap.set_bad('0.75') #75% grey
    vmin= -1.5
    vmax= 1.5
    #
    # tell the colormap what the maximum and minimum
    # values are
    #
    the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
    chlor=ma.array(chlor,mask=np.isnan(chlor))
    # create figure, add axes
    plt.close('all')
    fig=plt.figure(figsize=(12, 12))
    ax=fig.add_subplot(111)
    ## define parallels and meridians to draw.
    parallels=np.arange(-90, 90, 5)
    meridians=np.arange(0, 360, 5)
    proj.drawparallels(parallels, labels=[1, 0, 0, 0],\
                      fontsize=10, latmax=90)
    proj.drawmeridians(meridians, labels=[0, 0, 0, 1],\
                      fontsize=10, latmax=90)
    # draw coast & fill continents
    #map.fillcontinents(color=[0.25, 0.25, 0.25], lake_color=None) # coral
    out=proj.drawcoastlines(linewidth=1.5, linestyle='solid', color='k')
    x, y=proj(lon,lat)
    CS=proj.pcolormesh(x, y, np.log10(chlor), cmap=cmap,norm=the_norm)
    CBar=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
    CBar.set_label('log10 (mg/m^3) of the chlorophyll a concentration')
    title=ax.set_title('Chlorophyll-a {:s} to {:s} (grey=missing data)'.format(start,stop))
    fig.canvas.draw()
    plt.show()
