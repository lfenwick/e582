"""
   python translation of the perl script at
   http://oceancolor.gsfc.nasa.gov/cms/L3Bins.html

   usage:
            import tile_calc
            test18=tile_calc.tile_calc(18)
            print(test18.latlon2tile(49.2333,-123.25))
            should output  tile=343

"""


# The following functions are based on the pseudocode found in Appendix A of:
# 
# Campbell, J.W., J.M. Blaisdell, and M. Darzi, 1995:
# Level-3 SeaWiFS Data Products: Spatial and Temporal Binning Algorithms.
# NASA Tech. Memo. 104566, Vol. 32,
# S.B. Hooker, E.R. Firestone, and J.G. Acker, Eds.,
# NASA Goddard Space Flight Center, Greenbelt, Maryland

from __future__ import division
import numpy as np 

class tile_calc(object):

    def __init__(self,numrows):
        """
          initalize the class with the number of latitude rows in the tiling system.
          for Modis ocean color numrows=4320 for standard 4.64 km resolution

          initialize lattile (vector of length numrows containing center latitude of each row)
                     numtile (vector of length numrows containing number of tiles in each row)
                     basetile (vector of length numrows containing tile number at the starting point of each row
        """
        self.numrows=numrows
        basetile=[1]
        numtile=[]
        lattile=[]
        for row in range(numrows):
            lattile.append(((row + 0.5)*180.0/self.numrows) - 90.0)
            numtile.append(int(2*self.numrows*np.cos(lattile[row]*np.pi/180.0) + 0.5))
            if row > 0:
                basetile.append(basetile[row-1] + numtile[row-1])
        self.basetile=np.array(basetile)
        self.numtile=np.array(numtile)
        self.lattile=np.array(lattile)
        self.tottiles = basetile[numrows - 1] + numtile[numrows - 1] - 1

    def lat2row(self,lat):
        row=int((90. + lat)*self.numrows/180.)
        if row > self.numrows:
            row = self.numrows - 1
        return row

    def rowlon2tile(self,row,lon):
        lon = self.constrain_lon(lon)
        col = int((lon + 180.0)*self.numtile[row]/360.0)
        if col >= numtile[row]:
            col = numtile[row] - 1
        return self.basetile[row] + col

    def latlon2tile(self,lat,lon):
        lat = self.constrain_lat(lat)
        lon = self.constrain_lon(lon)
        row = self.lat2row(lat)
        col = int((lon + 180.0)*self.numtile[row]/360.0)
        if col >= self.numtile[row]:
            col = self.numtile[row] - 1

        return self.basetile[row] + col

    def tile2latlon(self,tile):
        row = self.numrows - 1
        if tile < 1:
            tile = 1
        while tile < self.basetile[row]:
            row-=1
        clat = self.lattile[row]
        clon = 360.0*(tile - self.basetile[row] + 0.5)/self.numtile[row] - 180.0
        return clat,clon

    def tile2bounds(self,tile):
        row=self.numrows -1
        if tile < 1:
          tile=1
        while tile < self.basetile[row]:
            row-=1
        north = self.lattile[row] + (90.0/self.numrows)
        south = self.lattile[row] - (90.0/self.numrows)
        lon = 360.0*(tile - self.basetile[row] + 0.5)/self.numtile[row] - 180.0
        west = lon - 180.0/self.numtile[row]
        east = lon + 180.0/self.numtile[row]
        return north,south,west,east

    def constrain_lat(self,lat):
        if lat > 90.:
            lat = 90
        if lat < -90.:
            lat = -90
        return lat

    def constrain_lon(self,lon):
        if lon < -180:
          lon += 360
        if lon > 180:
          lon -= 360
        return lon

if __name__ == "__main__":
    
    import tile_calc
    print('lat lon bounds for tile 367: ',test18.tile2bounds(367))
    print('bounds for 411: ',test18.tile2bounds(411))
    print('bounds for 412: ',test18.tile2bounds(412))
    print('central lat lon for 367: ',test18.tile2latlon(367))
    print('tile for UBC: ',test18.latlon2tile(49.2333,-123.25))
          
