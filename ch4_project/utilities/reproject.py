"""
   bin modis data into regular latitude and longitude bins   
"""

from __future__ import division
import numpy as np
import numba

def reproj_L1B(raw_data,bad_value,raw_lon, raw_lat, lonlim, latlim, res):
    
    '''
    =========================================================================================
    Reproject MODIS L1B file to a regular grid
    -----------------------------------------------------------------------------------------
    d_array, lon_array, lat_array, bin_count = reproj_L1B(raw_data,bad_value, raw_lon, raw_lat, lonlim, latlim, res)
    -----------------------------------------------------------------------------------------
    Input:
            raw_data: L1B data, N*M 2-D array.
            raw_lon: longitude info. N*M 2-D array.
            raw_lat: latitude info. N*M 2-D array.
            lonlim: range of longitude, a list.
            latlim: range of latitude, a list.
            res: resolution, single value.
    Output:
            d_array: L1B reprojected data.
            lon_array: reprojected longitude.
            lat_array: reprojected latitude.
            bin_count: how many raw data point included in a reprojected grid.
    Note:
            function does not perform well if "res" is larger than the resolution of input data.
            shape of "raw_data", "raw_lon", "raw_lat" must agree.
    =========================================================================================
    '''
    import numpy as np

    #
    # what do these lines do?
    #
    lon_edges=np.arange(lonlim[0], lonlim[1], res)
    lon_centers=np.empty_like(lon_edges)
    lon_centers[:-1]=(lon_edges[1:] + lon_edges[0:-1])/2.
    lon_centers[-1]=lon_centers[-2] + res
    lat_edges=np.arange(latlim[0], latlim[1], res)
    lat_centers=np.empty_like(lat_edges)
    lat_centers[:-1]=(lat_edges[1:] + lat_edges[0:-1])/2.
    lat_centers[-1]=lat_centers[-2] + res
    #print(np.diff(lon_centers[-5:]),np.diff(lat_centers[-5:]))
    #
    # meshgrid returns a lon_array shape=[len(lon_centers),len(lat_centers)],with lon_values down rows
    # and lat_array with lat values constant across columns.  These are the coordinates of
    #
    lon_array,lat_array=np.meshgrid(lon_centers,lat_centers)
    #print("here: ",lon_array,lat_array)
    #
    # for every value in raw_lon, find the index in lon_indices that
    # is the first index larger than that value (i.e. the right side
    # of the longitude bin
    #
    lon_indices=np.searchsorted(lon_edges, raw_lon.flat, 'right')
    lat_indices=np.searchsorted(lat_edges, raw_lat.flat, 'right')
    print('through search sorted')
    d_array=np.zeros([len(lat_edges), len(lon_edges)], dtype=np.float)
    bin_count=np.zeros([len(lat_edges), len(lon_edges)], dtype=np.int)
    
    for n in range(len(lat_indices)): #lat_indices and lon_indices both size of raw data
        bin_row=lat_indices[n]-1 # '-1' to get the index of the left edge of the bin
        bin_col=lon_indices[n]-1
        #
        # if the data is flagged as missing, assign that bin cell a value of np.nan
        # a single np.nan in the bin will cause that bin to be skipped subsequently
        #
        if (raw_data.flat[n] != bad_value) and (not np.isnan(d_array[bin_row, bin_col])):
            d_array[bin_row, bin_col] += raw_data.flat[n]
            bin_count[bin_row, bin_col] += 1
        elif np.isnan(d_array[bin_row, bin_col]):
            continue
        else:
            d_array[bin_row, bin_col] = np.nan
            bin_count[bin_row, bin_col]=0
            
    for i in range(lon_array.shape[0]):
        for j in range(lon_array.shape[1]):
            if bin_count[i, j] > 0:
                d_array[i, j]=d_array[i, j]/bin_count[i, j] 
            else:
                d_array[i, j]=np.nan
    return d_array, lon_array, lat_array, bin_count


def reproj_numba(raw_data,bad_value,raw_lon, raw_lat, lonlim, latlim, res):
    
    '''
    =========================================================================================
    Reproject MODIS L1B file to a regular grid
    -----------------------------------------------------------------------------------------
    d_array, lon_array, lat_array, bin_count = reproj_L1B(raw_data,bad_value, raw_lon, raw_lat, lonlim, latlim, res)
    -----------------------------------------------------------------------------------------
    Input:
            raw_data: L1B data, N*M 2-D array.
            raw_lon: longitude info. N*M 2-D array.
            raw_lat: latitude info. N*M 2-D array.
            lonlim: range of longitude, a list.
            latlim: range of latitude, a list.
            res: resolution, single value.
    Output:
            d_array: L1B reprojected data.
            lon_array: reprojected longitude.
            lat_array: reprojected latitude.
            bin_count: how many raw data point included in a reprojected grid.
    Note:
            function does not perform well if "res" is larger than the resolution of input data.
            shape of "raw_data", "raw_lon", "raw_lat" must agree.
    =========================================================================================
    '''
    import numpy as np

    #
    # what do these lines do?
    #
    lon_edges=np.arange(lonlim[0], lonlim[1], res)
    lon_centers=np.empty_like(lon_edges)
    lon_centers[:-1]=(lon_edges[1:] + lon_edges[0:-1])/2.
    lon_centers[-1]=lon_centers[-2] + res
    lat_edges=np.arange(latlim[0], latlim[1], res)
    lat_centers=np.empty_like(lat_edges)
    lat_centers[:-1]=(lat_edges[1:] + lat_edges[0:-1])/2.
    lat_centers[-1]=lat_centers[-2] + res
    #print(np.diff(lon_centers[-5:]),np.diff(lat_centers[-5:]))
    #
    # meshgrid returns a lon_array shape=[len(lon_centers),len(lat_centers)],with lon_values down rows
    # and lat_array with lat values constant across columns.  These are the coordinates of
    #
    lon_array,lat_array=np.meshgrid(lon_centers,lat_centers)
    #print("here: ",lon_array,lat_array)
    #
    # for every value in raw_lon, find the index in lon_indices that
    # is the first index larger than that value (i.e. the right side
    # of the longitude bin
    #
    lon_indices=np.searchsorted(lon_edges, raw_lon.flat, 'right')
    lat_indices=np.searchsorted(lat_edges, raw_lat.flat, 'right')
    print('through search sorted')
    d_array=np.zeros([len(lat_edges), len(lon_edges)], dtype=np.float)
    bin_count=np.zeros([len(lat_edges), len(lon_edges)], dtype=np.int)
    raw_lon_flat,raw_lat_flat=raw_lon.reshape(-1).astype(np.float32),raw_lat.reshape(-1).astype(np.float32)
    raw_data_flat=raw_data.reshape(-1).astype(np.float32)
    d_array, lon_array, lat_array, bin_count=fill_counts(lat_indices,lon_indices,
                  raw_data_flat,bad_value,
                  bin_count,d_array,lon_array,lat_array)

    return (d_array,lon_array,lat_array,bin_count)

@numba.jit(nopython=True)
def fill_counts(lat_indices,lon_indices,
                raw_data,bad_value,
                bin_count,d_array,lon_array,lat_array):

    for n in range(len(lat_indices)): #lat_indices and lon_indices both size of raw data
        bin_row=lat_indices[n]-1 # '-1' to get the index of the left edge of the bin
        bin_col=lon_indices[n]-1
        #
        # if the data is flagged as missing, assign that bin cell a value of np.nan
        # a single np.nan in the bin will cause that bin to be skipped subsequently
        #
        if (raw_data[n] != bad_value) and (not np.isnan(d_array[bin_row, bin_col])):
            d_array[bin_row, bin_col] += raw_data[n]
            bin_count[bin_row, bin_col] += 1
        elif np.isnan(d_array[bin_row, bin_col]):
            continue
        else:
            d_array[bin_row, bin_col] = np.nan
            bin_count[bin_row, bin_col]=0
            
    for i in range(lon_array.shape[0]):
        for j in range(lon_array.shape[1]):
            if bin_count[i, j] > 0:
                d_array[i, j]=d_array[i, j]/bin_count[i, j] 
            else:
                d_array[i, j]=np.nan
    return (d_array, lon_array, lat_array, bin_count)

