#!/usr/bin/env python
"""
  read the level1b file, geom file, cloud model file and the cloud mask file from dump_cloudmask.py and
  produce plots of the the particle phase

  usage:

  ./plot_palettes.py names.json

"""
from __future__ import division
import argparse
import h5py
import glob
from matplotlib import pyplot as plt
import site
site.addsitedir('../utilities')
from reproject import reproj_numba
import planck
import os,errno
import seaborn as sns
#
# compat module redefines importlib.reload if we're
# running python3
#
from compat import cpreload as reload
from planck import planckInvert
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
import textwrap
import io,json
from collections import OrderedDict as od
import compat
from compat import text_
import numpy.ma as ma
#
# two different ways of generating a colormap
# use this for continuous variables
from matplotlib.colors import LinearSegmentedColormap
#
# and use this for discrete colors
#
from matplotlib.colors import ListedColormap
import pdb


def make_plot(lcc_values):
    """
      set up the basic map projection details with coastlines and meridians
      return the projection object for further plotting
    """
    proj = Basemap(**lcc_values)
    parallels = np.arange(-90, 90, 1)
    meridians = np.arange(0, 360, 2)
    proj.drawparallels(parallels, labels=[1, 0, 0, 0],
                       fontsize=10, latmax=90)
    proj.drawmeridians(meridians, labels=[0, 0, 0, 1],
                       fontsize=10, latmax=90)
    # draw coast & fill continents
    # map.fillcontinents(color=[0.25, 0.25, 0.25], lake_color=None) # coral
    proj.drawcoastlines(linewidth=3., linestyle='solid', color='r')
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


def get_channels(l1b_file,geom_file,cloud_mask,mod06_file):
    with h5py.File(geom_file) as geom_h5,h5py.File(l1b_file) as l1b_h5:
        #channel31 is emissive channel 10
        #pdb.set_trace()
        band_names=text_(l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'].attrs['band_names'])
        band_names=band_names.split(',')
        index31=band_names.index('31')
        chan31=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'][index31,:,:]
        scale=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'].attrs['radiance_scales'][index31]
        offset=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'].attrs['radiance_offsets'][index31]
        chan31=(chan31 - offset)*scale
        index1=0  #channel 1 is first 250 meter reflective channel
        reflective=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_250_Aggr1km_RefSB'][0,:,:]
        scale=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_250_Aggr1km_RefSB'].attrs['reflectance_scales']
        offset=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_250_Aggr1km_RefSB'].attrs['reflectance_offsets']
        chan1=(reflective - offset[0])*scale[0]
        index29=band_names.index('29')

        chan29=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'][index29,:,:]
        scale=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'].attrs['radiance_scales'][index29]
        offset=l1b_h5['MODIS_SWATH_Type_L1B']['Data Fields']['EV_1KM_Emissive'].attrs['radiance_offsets'][index29]
        chan29=(chan29 - offset)*scale

        the_lon=geom_h5['MODIS_Swath_Type_GEO']['Geolocation Fields']['Longitude'][...]
        the_lat=geom_h5['MODIS_Swath_Type_GEO']['Geolocation Fields']['Latitude'][...]

    c31_bright=planckInvert(11.03,chan31)
    c29_bright=planckInvert(8.55,chan29)

    with h5py.File(cloud_mask) as cm_h5:
         maskout=cm_h5['cloudmask'][...]
         landout=cm_h5['landmask'][...]

    with h5py.File(mod06_file) as m06_h5:
        phase=m06_h5['mod06']['Data Fields']['Cloud_Phase_Infrared_1km'][...]
        phase=phase.astype(np.int32)
    return c31_bright,c29_bright,maskout,landout,phase,the_lon,the_lat

    
if __name__ == "__main__":

    linebreaks=argparse.RawTextHelpFormatter
    descrip=textwrap.dedent(globals()['__doc__'])
    parser = argparse.ArgumentParser(formatter_class=linebreaks,description=descrip)
    parser.add_argument('initfile',type=str,help='name of json file with filenames')
    args=parser.parse_args()

    with io.open(args.initfile,'r',encoding='utf8') as f:
        name_dict=json.loads(f.read(),object_pairs_hook=od)

    keys=['l1b_file', 'geom_file', 'mask_file', 'm06_file']

    the_files=[name_dict[key] for key in keys]

    plot_dir='plots'
    try:
        os.makedirs(plot_dir)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass  #not a problem if file exists

    c31_bright,c29_bright,maskout,landout,phase,the_lon,the_lat=\
      get_channels(*the_files)
    
    lcc_values,lon_res,lat_res=find_corners(the_lon,the_lat)
    lcc_values['fix_aspect']=True
    lcc_values['resolution']='c'
    lcc_values['projection']='lcc'


    plt.close('all')
    #
    # grid the data
    #
    missing_val=-999.  
    latlim=[lcc_values['llcrnrlat'],lcc_values['urcrnrlat']]
    lonlim=[lcc_values['llcrnrlon'],lcc_values['urcrnrlon']]
    res=0.027
    mask_grid, longrid, latgrid, bin_count = reproj_numba(maskout,missing_val, the_lon, the_lat, lonlim, latlim, res)
    land_grid, longrid, latgrid, bin_count = reproj_numba(landout,missing_val, the_lon, the_lat, lonlim, latlim, res)
    land_grid, longrid, latgrid, bin_count = reproj_numba(landout,missing_val, the_lon, the_lat, lonlim, latlim, res)
    c31bright_grid, longrid, latgrid, bin_count = reproj_numba(c31_bright,missing_val, the_lon, the_lat, lonlim, latlim, res)
    c29bright_grid, longrid, latgrid, bin_count = reproj_numba(c29_bright,missing_val, the_lon, the_lat, lonlim, latlim, res)
    #
    # need to change all 6's to 4's so the color map goes uniformly from 0-4
    #
    hit=phase==6
    phase[hit]=4
    phase_grid,longrid, latgrid, bin_count = reproj_numba(phase,missing_val, the_lon, the_lat, lonlim, latlim, res)
    #
    # plot the ungridded cloud model phase variable (MYD06_L2) in 5 discrete colors
    # using ListedColorMap 
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    #
    # make a 5 color palette
    #
    colors = ["royal blue", "baby blue", "eggshell", "burnt red", "soft pink"]
    print([the_color for the_color in colors])
    colors=[sns.xkcd_rgb[the_color] for the_color in colors]
    pal=ListedColormap(colors,N=5)
    #
    # the A2014127.2110 scene is a descending orbit, so south is on top
    # and west is on the right, need to rotate through 180 degrees
    #
    phase_rot=np.rot90(phase,2)
    CS=ax.imshow(phase_rot,cmap=pal)
    ax.set_title('ungridded phase map with 2 rotations')
    cax=fig.colorbar(CS)
    labels='0 -- cloud free,1 -- water cloud,2 -- ice cloud,3 -- mixed phase cloud,4 -- undetermined phase'
    labels=labels.split(',')
    ends=np.linspace(0,4,6)
    centers=(ends[1:] + ends[:-1])/2.
    cax.set_ticks(centers)
    cax.set_ticklabels(labels)
    fig.canvas.draw()
    fig.savefig('{}/ungridded_phase_map.png'.format(plot_dir))
    #
    # histogram the raw phase values
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    ax.hist(phase.ravel())
    ax.set_title('Raw: Gridded: Mask - 0 = Cloud,1 = 66% prob.\n Clear,2 = 95% prob. Clear,3 = 99% prob. Clear')
    fig.savefig('{}/ungridded_phase_hist.png'.format(plot_dir))
    #
    #  here is the gridded phase, lcc projection
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    lcc_values['ax']=ax
    proj=make_plot(lcc_values)
    x,y=proj(longrid,latgrid)
    pal.set_over('r')
    pal.set_under('k')
    pal.set_bad('0.75') #75% grey
    phase_grid=ma.array(phase_grid,mask=np.isnan(phase_grid))
    CS=proj.ax.pcolormesh(x,y,phase_grid,cmap=pal)
    cax=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
    labels='0 -- cloud free,1 -- water cloud,2 -- ice cloud,3 -- mixed phase cloud,4 -- undetermined phase'
    labels=labels.split(',')
    #
    # label the colorbar axis with centered values
    #
    ends=np.linspace(0,4,6)
    centers=(ends[1:] + ends[:-1])/2.
    cax.set_ticks(centers)
    cax.set_ticklabels(labels)
    cax.set_label('phase')
    proj.ax.set_title('phase')
    proj.ax.figure.canvas.draw()
    fig.savefig('{}/gridded_phase_map.png'.format(plot_dir))
    #
    # 
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    ax.hist(phase_grid.compressed())
    ax.set_title('Gridded: Mask - 0 = Cloud,1 = 66% prob.\n Clear,2 = 95% prob. Clear,3 = 99% prob. Clear')
    fig.savefig('{}/gridded_phase_hist.png'.format(plot_dir))
    #
    # compare this to the 8 - 11 micron brightness temperaure difference
    # this is a continuous variable so use 256 colors, 
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    colors=sns.color_palette('coolwarm')
    pal=LinearSegmentedColormap.from_list('test',colors)
    pal.set_bad('0.75') #75% grey
    pal.set_over('r')
    pal.set_under('k')
    vmin= -5.
    vmax= 5.
    the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
    lcc_values['ax']=ax
    proj=make_plot(lcc_values)
    tdiff_grid=c29bright_grid - c31bright_grid
    tdiff_grid=ma.array(tdiff_grid,mask=np.isnan(tdiff_grid))
    CS=proj.ax.pcolormesh(x,y,tdiff_grid,cmap=pal,norm=the_norm)
    cax=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
    ax.set_title('TB 8 miron - TB 11 micron')
    fig.savefig('{}/gridded_TBdiff_map.png'.format(plot_dir))
    #
    # now plot the cloud mask from MYD35_L2 for compairson
    #
    fig,ax=plt.subplots(1,1,figsize=(12,12))
    colors = ["eggshell","baby blue","sky blue","royal blue"]
    colors=[sns.xkcd_rgb[the_color] for the_color in colors]
    colors=sns.color_palette("Blues")
    colors=sns.color_palette('coolwarm_r')
    pal=ListedColormap(colors,N=4)
    pal.set_bad('0.75') #75% grey
    labels='0:Cloud,1:66% Clear,2:95% Clear,3: 99% Clear'
    labels=labels.split(',')
    ends=np.linspace(0,3,5)
    centers=(ends[1:] + ends[:-1])/2.
    lcc_values['ax']=ax
    proj=make_plot(lcc_values)
    mask_grid=ma.array(mask_grid,mask=np.isnan(mask_grid))
    CS=proj.ax.pcolormesh(x,y,mask_grid,cmap=pal)
    cax=proj.colorbar(CS, 'right', size='5%', pad='5%')
    cax.set_ticks(centers)
    cax.set_ticklabels(labels)
    proj.ax.set_title('Mask - 0 = Cloud,1 = 66% prob. Clear,2 = 95% prob. Clear,3 = 99% prob. Clear')
    proj.ax.figure.canvas.draw()
    fig.savefig('{}/mapped_cloudmask.png'.format(plot_dir))
    #
    # finally, for OSX and linux, read the high cloud bit that uses the 6.7 micron
    # brightness temperature
    #
    try:
        import bitmap
        rawmask_h5=h5py.File(name_dict['rawmask_file'])
        all_mask=rawmask_h5['mod35']['Data Fields']['Cloud_Mask'][...]
        #get bit 7 of the second byte
        #http://modis-atmos.gsfc.nasa.gov/_specs/MOD35_L2.CDL.fs
        high_cloud_flag=bitmap.getbit(all_mask[1,:,:],7)
        
        fig,ax=plt.subplots(1,1,figsize=(12,12))
        ax.hist(high_cloud_flag.ravel())
        ax.set_title('histogrammed 6.7 micron raw bit')
        
        high_grid, longrid, latgrid, bin_count = reproj_numba(high_cloud_flag,missing_val, the_lon, the_lat, lonlim, latlim, res)
        high_grid=ma.array(high_grid,mask=np.isnan(high_grid))
        
        fig,ax=plt.subplots(1,1,figsize=(12,12))
        ax.hist(high_grid.compressed())
        ax.set_title('histogram gridded 6.7 micron  bit')
        fig.savefig('{}/hist_gridded_bit.png'.format(plot_dir))
        #
        # now map this using a two color palette
        #
        colors = ["eggshell","sky blue"]
        colors=[sns.xkcd_rgb[the_color] for the_color in colors]
        pal=ListedColormap(colors,N=2)
        pal.set_bad('0.75') #75% grey
        pal.set_over('r')
        pal.set_under('k')
        vmin= 0.
        vmax= 1.
        the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
        ends=np.linspace(0,1,3)
        centers=(ends[1:] + ends[:-1])/2.
        fig,ax=plt.subplots(1,1,figsize=(12,12))
        lcc_values['ax']=ax
        proj=make_plot(lcc_values)
        CS=proj.ax.pcolormesh(x,y,high_grid,cmap=pal,norm=the_norm)
        cax=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
        labels='0:yes,1:no'
        labels=labels.split(',')
        cax.set_ticks(centers)
        cax.set_ticklabels(labels)
        ax.set_title('6.7 micron high cloud flag')
        fig.savefig('{}/map_gridded_bit.png'.format(plot_dir))
    except ImportError:
        print("couldn't import bitmap, skipping high cloud mask plot")
   
    
    plt.show()

