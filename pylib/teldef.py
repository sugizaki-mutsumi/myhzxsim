#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits

from scipy.spatial.transform import Rotation as R


class Teldef() :

    def __init__(self, caldbfname=None):


        fptr = pyfits.open(caldbfname)
        hdr = fptr[0].header

        self.filename = hdr['FILENAME']

        self.telescop = hdr['TELESCOP']
        self.instrume = hdr['INSTRUME']
        self.detnam   = hdr['DETNAM']

        self.raw_xsiz = hdr['RAW_XSIZ']  ### num of pixels
        self.rawxpix1 = hdr['RAWXPIX1']  ### 1st pixel (0 or 1)
        self.raw_xscl = hdr['RAW_XSCL']  ### pixel size
        self.raw_xcol = hdr['RAW_XCOL']  ### column name
        self.raw_ysiz = hdr['RAW_YSIZ']  ### num of pixels
        self.rawypix1 = hdr['RAWYPIX1']  ### 1st pixel (0 or 1)
        self.raw_yscl = hdr['RAW_YSCL']  ### pixel size
        self.raw_ycol = hdr['RAW_YCOL']  ### column name
        self.raw_unit = hdr['RAW_UNIT']  ### units        

        self.det_xsiz = hdr['DET_XSIZ']  ### num of pixels
        self.detxpix1 = hdr['DETXPIX1']  ### 1st pixel (0 or 1)
        self.det_xscl = hdr['DET_XSCL']  ### pixel size   
        self.det_xcol = hdr['DET_XCOL']  ### column name  
        self.det_ysiz = hdr['DET_YSIZ']  ### num of pixels
        self.detypix1 = hdr['DETYPIX1']  ### 1st pixel (0 or 1)
        self.det_yscl = hdr['DET_YSCL']  ### pixel size   
        self.det_ycol = hdr['DET_YCOL']  ### column name  
        self.det_unit = hdr['DET_UNIT']  ### units        

        self.sky_xsiz = hdr['SKY_XSIZ']  ### num of pixels
        self.skyxpix1 = hdr['SKYXPIX1']  ### 1st pixel    
        self.sky_xcol = hdr['SKY_XCOL']  ### column name  
        self.sky_ysiz = hdr['SKY_YSIZ']  ### num of pixels
        self.skyypix1 = hdr['SKYYPIX1']  ### 1st pixel    
        self.sky_ycol = hdr['SKY_YCOL']  ### column name  

        self.sky_unit = hdr['SKY_UNIT']  ### units
        #self.sky_from = hdr['SKY_FROM']  

        self.rawxpix_min = self.rawxpix1 -0.5                  ### float
        self.rawxpix_max = self.rawxpix1 + self.raw_xsiz -0.5  ### float
        self.rawypix_min = self.rawypix1 -0.5                  ### float
        self.rawypix_max = self.rawypix1 + self.raw_ysiz -0.5  ### float

        self.rawxpix_cen = self.rawxpix1 + 0.5*self.raw_xsiz -0.5  
        self.rawypix_cen = self.rawypix1 + 0.5*self.raw_ysiz -0.5  
        
        
        self.detxpix_min = self.detxpix1 -0.5                  ### float
        self.detxpix_max = self.detxpix1 + self.det_xsiz -0.5  ### float
        self.detypix_min = self.detypix1 -0.5                  ### float
        self.detypix_max = self.detypix1 + self.det_ysiz -0.5  ### float

        self.detxpix_cen = self.detxpix1 + 0.5*self.det_xsiz -0.5  
        self.detypix_cen = self.detypix1 + 0.5*self.det_ysiz -0.5  
        
        
        self.skyxpix_min = self.skyxpix1 -0.5                  ### float
        self.skyxpix_max = self.skyxpix1 + self.sky_xsiz -0.5  ### float
        self.skyypix_min = self.skyypix1 -0.5                  ### float
        self.skyypix_max = self.skyypix1 + self.sky_ysiz -0.5  ### float

        self.skyxpix_cen = self.skyxpix1 + 0.5*self.sky_xsiz -0.5  
        self.skyypix_cen = self.skyypix1 + 0.5*self.sky_ysiz -0.5  
        

        
        #self.alignm11 = hdr['ALIGNM11']
        #self.alignm12 = hdr['ALIGNM12']
        #self.alignm13 = hdr['ALIGNM13']
        #self.alignm21 = hdr['ALIGNM21']
        #self.alignm22 = hdr['ALIGNM22']
        #self.alignm23 = hdr['ALIGNM23']
        #self.alignm31 = hdr['ALIGNM31']
        #self.alignm32 = hdr['ALIGNM32']
        #self.alignm33 = hdr['ALIGNM33']

        m11 = hdr['ALIGNM11']
        m12 = hdr['ALIGNM12']
        m13 = hdr['ALIGNM13']
        m21 = hdr['ALIGNM21']
        m22 = hdr['ALIGNM22']
        m23 = hdr['ALIGNM23']
        m31 = hdr['ALIGNM31']
        m32 = hdr['ALIGNM32']
        m33 = hdr['ALIGNM33']

        #self.alignmat = np.mat([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]) 
        self.alignmat = np.asmatrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]) 
        self.alignrot = R.from_matrix(self.alignmat)  #.inv()

        self.rollsign = hdr['ROLLSIGN']
        self.rolloff  = hdr['ROLLOFF']

        
        self.focallen = hdr['FOCALLEN']
        self.optaxisx = hdr['OPTAXISX']
        self.optaxisy = hdr['OPTAXISY']

        if self.detnam[0]=="X" :
            self.detid  = int(self.detnam[2:])   ### Num of 1...16
        else : 
            self.detid  = int(self.detnam[4:])   ### Num of 1...48
        
if __name__=="__main__" :


    fname = sys.argv[1]
    teldef = Teldef(fname)

    if teldef!=None :
        
        print ("filename = %s"%(teldef.filename))
        print ("telescop = %s"%(teldef.telescop))
        print ("instrume = %s"%(teldef.instrume))
        print ("detnam   = %s"%(teldef.detnam))

        print ("detid = %d"%(teldef.detid))
        #print ("subid = %d"%(teldef.subid))

        print ("raw_xsiz = %d "%(teldef.raw_xsiz))
        print ("raw_xscl = %lf"%(teldef.raw_xscl))
        print ("rawxpix1 = %lf"%(teldef.rawxpix1))
        print ("alignmat  = \n", teldef.alignmat)
               
        
