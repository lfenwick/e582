#!/usr/bin/env python

from __future__ import division
import argparse
import h5py
import glob
from matplotlib import pyplot as plt
import site
site.addsitedir('./utilities')
site.addsitedir('./AIRS')
from reproject import reproj_numba
import planck
import os,errno
#
# compat module redefines importlib.reload if we're
# running python3

#
from compat import cpreload as reload
reload(planck)
from planck import planckInvert
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
import textwrap
import io,json
from collections import OrderedDict as od
import compat
reload(compat)
from compat import text_
import numpy.ma as ma

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


path = '/Users/lindsayfenwick/Google Drive/Active/my_notebooks/AIRS/*.h5'
keep_list=[]
for fname in glob.glob(path):
	with h5py.File(fname) as airs_h5:
    		missing_val = airs_h5['/L2_Standard_atmospheric&surface_product']['Data Fields']['CH4_total_column'].attrs['_FillValue']
    		ch4_col = airs_h5['/L2_Standard_atmospheric&surface_product']['Data Fields']['CH4_total_column'][...]
    		qc_ch4 = airs_h5['/L2_Standard_atmospheric&surface_product']['Data Fields']['CH4_total_column_QC'][...]
    		the_lon = airs_h5['/L2_Standard_atmospheric&surface_product']['Geolocation Fields']['Longitude'][...]
    		the_lat = airs_h5['/L2_Standard_atmospheric&surface_product']['Geolocation Fields']['Latitude'][...]
       # mixing_ratio = airs_h5['/L2_Standard_atmospheric&surface_product']['Data Fields']['CH4VMRSurf'][...]
    	bad_values=ch4_col == missing_val
    	ch4_col_converted = ch4_col*(16*1.66054e-27/1.e-4) #convert from molecules ch4/cm^2 to kg/m^2
        avg_ch4=np.median(ch4_col_converted)
	ch4_col_converted[bad_values]=missing_val
	keep_dict=dict(fname=fname,missing_val=missing_val,ch4_col=ch4_col_converted,the_lon=the_lon,the_lat=the_lat,
			qc_ch4=qc_ch4,bad_values=bad_values,avg_ch4=avg_ch4)
       	keep_list.append(keep_dict)

the_lon=keep_list[25]['the_lon']
the_lat=keep_list[25]['the_lat']
missing_val = keep_list[25]['missing_val']

lcc_values,lon_res,lat_res=find_corners(the_lon,the_lat)
lcc_values['fix_aspect']=True

lcc_values['resolution']='c'
lcc_values['projection']='lcc'

latlim=[lcc_values['llcrnrlat'],lcc_values['urcrnrlat']]
lonlim=[lcc_values['llcrnrlon'],lcc_values['urcrnrlon']]
res=0.5
ch4_col_grid, longrid, latgrid, bin_count = reproj_numba(keep_list[25]['ch4_col'],missing_val, the_lon, the_lat, lonlim, latlim, res)
#phase_grid,longrid, latgrid, bin_count = reproj_numba(phase,missing_val, the_lon, the_lat, lonlim, latlim, res)

    
cmap=cm.YlGn  #see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
cmap.set_over('r')
cmap.set_under('b') 
cmap.set_bad('0.75') #75% grey
vmin= 0.009
vmax= 0.011

the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
fig,ax=plt.subplots(1,1,figsize=(10,10))
    #
    # tell Basemap what axis to plot into
    #
lcc_values['ax']=ax
proj=make_plot(lcc_values)
x,y=proj(longrid,latgrid)
CS=proj.ax.pcolormesh(x,y,ch4_col_grid,cmap=cmap,norm=the_norm)
CBar=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
CBar.set_label('Methane Total Column')
proj.ax.set_title('Methane total column')
proj.ax.figure.canvas.draw()
plt.show()   

# for key in keep_list:
#   if key == 24: 
#     print key['ch4_col']

# print keep_list[26]['ch4_col']


# count_good_vals = 0
# for key in keep_list:
  
#   for array in key['qc_ch4']:
#     for number in array:
#       if number >= 1:
#         count_good_vals = count_good_vals + 1
# print count_good_vals

# average = []
# quality = []
# guess = []
# for key in keep_list:
#   average.append(key['avg_ch4']) #average ch4 for each kernel
#   bad_vals = key['qc_ch4'] == 0
#   quality.append(key['qc_ch4'])
#   guess.append(key['first_guess'])

# plt.hist(guess)

#plt.hist(quality, bins = []
#plt.show()

#plt.hist(keep_dict['qc_ch4'])
#plt.show()


#plt.hist(all_ch4)
#plt.show()

#I want to 

# # for keep_dict in keep_list:
#   plt.hist(keep_dict[qc_ch4])
# plt.show()


#plt.scatter(the_lon,the_lat)
#plt.hist(keep_list('ch4_col_converted')
#plt.show()


#fig,ax=plt.subplots(1,1,figsize=(10,10))
#ax.hist(ch4_col_converted.ravel())
#ax.set_title('pre-gridded')

#plt.scatter(the_lon,the_lat)
#plt.show()

# lcc_values,lon_res,lat_res=find_corners(the_lon,the_lat)
# lcc_values['fix_aspect']=True
# ##lcc_values['lat_0']=55.
# ##lcc_values['lat_1']=90.
# ##lcc_values['lat_2']=55.
# ##lcc_values['llcrnrlat']=55.
# ##lcc_values['urcrnrlat']=60.
# ##lcc_values['llcrnrlon']= 0.
# ##lcc_values['urcrnrlon']= 360.
# ##lcc_values['width']=200.e3
# ##lcc_values['height']=500.e3
# lcc_values['resolution']='c'
# lcc_values['projection']='lcc'

# latlim=[lcc_values['llcrnrlat'],lcc_values['urcrnrlat']]
# lonlim=[lcc_values['llcrnrlon'],lcc_values['urcrnrlon']]
# res=0.5
# ch4_col_grid, longrid, latgrid, bin_count = reproj_numba(ch4_col_converted,missing_val, the_lon, the_lat, lonlim, latlim, res)
# #phase_grid,longrid, latgrid, bin_count = reproj_numba(phase,missing_val, the_lon, the_lat, lonlim, latlim, res)

    
# cmap=cm.YlGn  #see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
# cmap.set_over('r')
# cmap.set_under('b') 
# cmap.set_bad('0.75') #75% grey
# vmin= 0.009
# vmax= 0.011

# the_norm=Normalize(vmin=vmin,vmax=vmax,clip=False)
# fig,ax=plt.subplots(1,1,figsize=(10,10))
#     #
#     # tell Basemap what axis to plot into
#     #
# lcc_values['ax']=ax
# proj=make_plot(lcc_values)
# x,y=proj(longrid,latgrid)
# CS=proj.ax.pcolormesh(x,y,ch4_col_grid,cmap=cmap,norm=the_norm)
# CBar=proj.colorbar(CS, 'right', size='5%', pad='5%',extend='both')
# CBar.set_label('Methane Total Column')
# proj.ax.set_title('Methane total column')
# proj.ax.figure.canvas.draw()
# plt.show()   