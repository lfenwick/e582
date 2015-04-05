#cython: embedsignature=True
import numpy as np
cimport numpy as np
from libc.stdint cimport int8_t

#try from libc.stdint cimport int64_t 

cdef extern  void readcloud_cpp(int8_t* byteone,int8_t* maskout,int nvals)

cdef extern  void readland_cpp(int8_t* byteone,int8_t* landout,int nvals)

cdef extern  void readthin_cirrus_cpp(int8_t* byteone,int8_t* highout,int nvals)

cdef extern  void readhigh_cloud_cpp(int8_t* byteone,int8_t* thinout,int nvals)

cdef extern  void readbit(int8_t* inbyte,int8_t* outbyte,int bitnum, int nvals)

def getmask_zero(object byteone):
    """
       http://modis-atmos.gsfc.nasa.gov/MOD35_L2/format.html
       http://modis-atmos.gsfc.nasa.gov/MOD35_L2/index.html
       http://modis-atmos.gsfc.nasa.gov/_specs/MOD35_L2.CDL.fs
       input is 2 dimensional numpy array of type np.int8
       containing 0 byte of the cloud mask

       output is a tuple of two 2 dimensional np.int8 arrays

       maskout has each pixel's cloud probability:

        0 = Cloud            
        1 = 66% prob. Clear  
        2 = 95% prob. Clear  
        3 = 99% prob. Clear

       landout has land/water values

       0=Water   
       1=Coastal 
       2=Desert  
       3=Land    

    """
    byteone=np.ascontiguousarray(byteone)
    cdef int nvals= byteone.size
    saveShape=byteone.shape
    cdef np.ndarray[np.int8_t,ndim=2] c_byte=byteone
    cdef int8_t* dataPtr=<int8_t*> c_byte.data
    cdef np.ndarray[np.int8_t,ndim=2] maskout
    maskout=np.empty(saveShape,dtype=np.int8)
    cdef int8_t* maskPtr=<int8_t*> maskout.data
    cdef np.ndarray[np.int8_t,ndim=2] landout=np.empty(saveShape,dtype=np.int8)
    cdef int8_t* landPtr=<int8_t*> landout.data
    readcloud_cpp(dataPtr, maskPtr,nvals)
    readland_cpp(dataPtr, landPtr,nvals)
    out=(maskout,landout)
    return out


def getmask_one(object byteone):
    """
       input is 2 dimensional numpy array of type np.int8
       containing  byte 1 of the cloud mask

       output is a tuple of two 2 dimensional np.int8 arrays

       thinout has 0 if thin cirrus, 1 if not, so clear=True


       highout has 0 if high cloud, 1 if not, so clear=True

    """
    byteone=np.ascontiguousarray(byteone)
    cdef int nvals= byteone.size
    saveShape=byteone.shape
    cdef np.ndarray[np.int8_t,ndim=2] c_byte=byteone
    cdef np.int8_t* dataPtr=<np.int8_t*> c_byte.data
    cdef np.ndarray[np.int8_t,ndim=2] thinout
    thinout=np.empty(saveShape,dtype=np.int8)
    cdef np.int8_t* thinPtr=<np.int8_t*> thinout.data
    cdef np.ndarray[np.int8_t,ndim=2] highout=np.empty(saveShape,dtype=np.int8)
    cdef np.int8_t* highPtr=<np.int8_t*> highout.data
    readthin_cirrus_cpp(dataPtr, thinPtr,nvals)
    readhigh_cloud_cpp(dataPtr, highPtr,nvals)
    out=(thinout,highout)
    return out

def getbit(object inbyte,int bitnum):
    """
       inbyte is 2 dimensional numpy array of type np.int8
       containing  one of the six bytes from /mod35/Data Fields/Cloud_Mask

       bitnum is the bit number (0-7) to read out of that byte

       returns:

           an int8 numarray of 0,1 values with the same shape
           as inbyte

    """
    pybute=np.ascontiguousarray(inbyte)
    cdef int nvals= inbyte.size
    cdef np.int8_t [:,:] inbyte_view=inbyte
    cdef np.int8_t* dataPtr=<np.int8_t*> &inbyte_view[0,0]
    bitout=np.empty_like(inbyte)
    cdef np.int8_t [:,:] bitout_view=bitout
    cdef np.int8_t* bitPtr= <np.int8_t*> &bitout_view[0,0]
    readbit(dataPtr, bitPtr,bitnum,nvals)
    return bitout



