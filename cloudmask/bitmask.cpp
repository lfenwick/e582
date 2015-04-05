//#include <tr1/memory>
#include <bitset>
#include <iostream>

/*
http://modis-atmos.gsfc.nasa.gov/MOD35_L2/format.html
byte 0
+------+------------------------------+------------------------+
| bits |       Field Description      | Bit Interpretation Key |
+------+------------------------------+------------------------+
| 0    | Cloud Mask Flag              | 0 = Not Determined     |
|      |                              | 1 = Determined         |
| 1-2  | Unobstructed FOV Quality Flag| 0 = Cloud              | 
|      |                              | 1 = 66% prob. Clear    |
|      |                              | 2 = 95% prob. Clear    |
|      |                              | 3 = 99% prob. Clear    | 
| 3    | Day/Night Flag               | 0 = Night              |
|      |                              | 1 = Day                | 
| 4    | Sunglint Flag                | 0 = Yes                |
|      |                              | 1 = No                 |
| 5    | Snow/Ice Background Flag     | 0 = Yes                |
|      |                              | 1 = No                 |
| 6-7  | Land/Water Background Flag   | 0=Water                |
|      |                              | 1=Coastal              |
|      |                              | 2=Desert               |
|      |                              | 3=Land                 |
+------+------------------------------+------------------------+

byte 1

0        Non-Cloud Obstruction Flag       0 = Yes 
                                          1 = No

1       Thin Cirrus Detected (Solar Test) 0 = Yes
                                          1 = No

2       Shadow Flag                       0 = Yes
                                          1 = No

3       Thin Cirrus Detected (Infrared Test)    0 = Yes
                                                1 = No
4       Adjacent Cloud Detected
(within surrounding 1 km pixels)        0 = Yes
                                        1 = No

5       Cloud Flag (IR Threshold Test)  0 = Yes
                                        1 = No

6       High Cloud Flag (CO2 Test)      0 = Yes
                                        1 = No

7       High Cloud Flag (6.7 micron Test) 0 = Yes
                                          1 = No 
byte 3

0        High Cloud Flag (1.38 micron Test)      0 = Yes
                                                 1 = No

1       High Cloud Flag (3.7-12 micron Test)    0 = Yes
                                                1 = No

2       Cloud Flag (IR Temperature Difference Test)     0 = Yes
                                                        1 = No

3       Cloud Flag (3.7-11 micron Test)         0 = Yes
                                                1 = No

4       Cloud Flag (Visible Reflectance Test)   0 = Yes
                                                1 = No

5       Cloud Flag (Visible Reflectance Ratio Test)     0 = Yes
                                                        1 = No

6       Cloud Flag (.935/.87 Reflectance Test)  0 = Yes
                                                1 = No

7       Cloud Flag (3.7-3.9 micron Test)        0 = Yes
                                                1 = No 



*/


extern "C" {
  void readcloud_cpp(char* byteone,char* maskout,int nvals){
    //byte 0, bits 1-2
    for(int i=0; i < nvals; ++i){
      std::bitset<8> myseq(byteone[i]);
      maskout[i] = myseq[1] +(myseq[2]*2);
    }
  }
  void readland_cpp(char* byteone,char* landout,int nvals){
    //byte 0, bit 6-7
    for(int i=0; i < nvals; ++i){
      std::bitset<8> myseq(byteone[i]);
      landout[i] = myseq[6]+(myseq[7]*2);
    }
  }
  void readthin_cirrus_cpp(char* byteone,char* cirrusout,int nvals){
    //byte 1, bit 3
    for(int i=0; i < nvals; ++i){
      std::bitset<8> myseq(byteone[i]);
      cirrusout[i] = myseq[3];
    }
  }
  void readhigh_cloud_cpp(char* byteone,char* highout,int nvals){
    //byte 1, bit 6
    for(int i=0; i < nvals; ++i){
      std::bitset<8> myseq(byteone[i]);
      highout[i] = myseq[6];
    }
  }

  void readbit(char* dataPtr,char* bitPtr,int bitnum, int nvals){
    for(int i=0; i < nvals; ++i){
      std::bitset<8> myseq(dataPtr[i]);
      bitPtr[i] = myseq[bitnum];
    }
  }

}





