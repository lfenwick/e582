# coding: utf-8
"""
   if write=True:
       read a NASA ocean color binned data file in hdf5 format
       calculate the mean chlorophyll concentration for every
       bin and create a pandas dataframe indexed on bin.  Save
       this dataframe for later use as a new hdf5 with the filename
       stored in the variable out_h5 ('store.h5')
   if write=False:
       read the dataframe in from out_h5
   then:
       calculate an array of lats and lons and find the mean
       chlorophyll concentration for each lat/lon pair
       write the lat,lon and chlor arrays out to a new
       hdf file with name stored in final_file ('navigated.h5')
"""

from __future__ import print_function, division
from builtins import range

import glob  #this module gets file names using wildcards
import site
site.addsitedir('../utilities')
import tile_calc
import h5py
import matplotlib.pyplot as plt
import numpy as np
import tile_calc
from contexttimer import Timer
import pickle
import pandas as pd
import time
import datetime as dt

print(__name__)
if __name__ == "__main__":
    #
    # store a dataframe with bins,meanvals for later use
    # in file out_h5
    # if write=True, or read from it if write=Fales
    #
    write=True
    out_h5='store.h5'
    #
    # this output file holds the gridded image
    #
    final_file='navigated.h5'
    nlats=500
    nlons=500
    lats=np.linspace(40.,50.,nlats)
    lons=np.linspace(-121,-131.,nlons)
    binned_file=glob.glob('../dataset/A20101522010181.L3b_MO_CHL.h5')[0]
    with  h5py.File(binned_file,'r') as infile:
        root_key=list(infile.keys())[0]
        #
        # turn day of year into a month and day
        # and save so we can write out as attributes
        # of our output files
        #
        start_day=int(infile.attrs['Start Day_GLOSDS'])  #convert from 16 bit to 64 bit int
        start_year=infile.attrs['Start Year_GLOSDS']
        start=dt.datetime(start_year-1,12,31) + dt.timedelta(days=start_day)
        end_day=int(infile.attrs['End Day_GLOSDS'])
        end_year=infile.attrs['End Year_GLOSDS']
        end=dt.datetime(end_year-1,12,31) + dt.timedelta(days=end_day)
        start_date=start.strftime('%Y-%m-%d')
        end_date=end.strftime('%Y-%m-%d')
        binlist=infile[root_key]['BinList']
        chlor_a=infile[root_key]['chlor_a']
        veclength=binlist.shape[0]
        print('reading length: ',veclength)
        chlor_a_data=chlor_a['chlor_a_sum'][:veclength]
        #chlor_a_sq_data=chlor_a['chlor_a_sum_sq'][:veclength]
        weights_data=binlist['weights'][:veclength]
        binnums=binlist['bin_num'][:veclength]
        out = np.empty((veclength,),dtype=[('binnum','>i4'),('chlor_a_mean','>f4')])
        #
        # keep track of elapsed time with the Timer object
        #
        with Timer() as t:
            if write:
                #
                # fill the structured array with tile,chlorophyll pairs
                #
                for i in range(veclength):
                    meanval=chlor_a_data[i]/weights_data[i]
                    out[i]=(binnums[i],chlor_a_data[i]/weights_data[i])
                print("time to create structured array: ",t.elapsed)
                #
                # create a pandas dataframe using the structured array
                # indexed by tile number
                #
                the_df=pd.DataFrame.from_records(out,index='binnum')
                print("time to create dataframe: ",t.elapsed)
                with pd.HDFStore('store.h5','w') as store:
                    store.put('chlor_a_mean',the_df,format='table')
                #
                # open the file a second time to write the attributes
                #
                with  h5py.File(out_h5,'a') as f:
                    f.attrs['history']='created by satelliteD.py'
                    f.attrs['created_on']=time.strftime("%c")
                    f.attrs['start_date']=start_date
                    f.attrs['end_date']=end_date
                    units='micrograms/m^3'
                    title='mean chlorophyll concentration'
                    f['/chlor_a_mean'].attrs['units']=units
                    f['/chlor_a_mean'].attrs['title']=title
                print('time to write dataframe: ',t.elapsed)
            else:
                #
                # reuse data from a previous file
                #
                with pd.HDFStore(out_h5,'r') as store:
                    the_df=store['chlor_a_mean']
                with  h5py.File(out_h5,'a') as f:
                    history=f.attrs['history']
                    created_on=f.attrs['created_on']
                    start_date=f.attrs['start_date']
                    end_date=f.attrs['end_date']
                    units=f['/chlor_a_mean'].attrs['units']
                    title=f['/chlor_a_mean'].attrs['title']

                    
        num_rows=4320
        tc=tile_calc.tile_calc(num_rows)
        tile_find=tile_calc.tile_calc(num_rows)
        lat_array=np.empty([nlats,nlons],dtype=np.float32)
        lon_array=np.empty_like(lat_array)
        chlor_array=np.empty_like(lat_array)
        with Timer() as t:
            with pd.HDFStore(out_h5,'r') as store:
                for row,the_lat in enumerate(lats):
                    for col,the_lon in enumerate(lons):
                        lat_array[row,col]=the_lat
                        lon_array[row,col]=the_lon
                        tile=tc.latlon2tile(the_lat,the_lon)
                        if tile in the_df.index:
                            #
                            # look up the chloropyll concentration
                            # using the bin index of the dataframe
                            #
                            chlor_array[row,col]=the_df.loc[tile]
                        else:
                            chlor_array[row,col]=np.nan
            print('through calc: ',t.elapsed)

            #
            # write the navigated data out to a new hdf file
            #
            with h5py.File(final_file, "w") as f:
                dset_lat = f.create_dataset("lattitude", lat_array.shape, dtype=lat_array.dtype)
                dset_lat[...]=lat_array[...]
                dset_lon = f.create_dataset("longitude", lon_array.shape, dtype=lon_array.dtype)
                dset_lon[...]=lon_array[...]
                dset_var = f.create_dataset("chlor_a_mean", chlor_array.shape, dtype=chlor_array.dtype)
                dset_var[...]=chlor_array[...]
                dset_var.attrs['units']=units
                dset_var.attrs['title']=title
                f.attrs['created_on']=time.strftime("%c")
                print('through final write: ',t.elapsed)
                f.attrs['history']='created by satelliteD.py'
                f.attrs['created_on']=time.strftime("%c")
                f.attrs['start_date']=start_date
                f.attrs['end_date']=end_date


        
            
            

                    

                

                
                


                
