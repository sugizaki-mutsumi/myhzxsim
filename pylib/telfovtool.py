#!/usr/bin/env python

import os, sys
import math
import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as wcs
from astropy.coordinates import SkyCoord

#from pyquaternion import Quaternion
from scipy.spatial.transform import Rotation as R

from teldef import Teldef
import telcoord

def gen_telfovreg(regfname, attrot, teldef_list, coordsys='icrs',
                  color_list=["white", "red", "green", "blue"]) :
    
    if os.path.isfile(regfname) :
        os.remove(regfname)
        
    regfile = open(regfname, "w")
    regfile.write('# Region file format: DS9 version 4.1\n')
    regfile.write('%s\n'%(coordsys))
        
    ngrid = 10
    idx=0

    ncolor = len(color_list)
    for teldef in teldef_list :
        colorname = color_list[idx%ncolor]

        xmin = teldef.detxpix1
        xmax = teldef.det_xsiz
        ymin = teldef.detypix1
        ymax = teldef.det_ysiz


        ### text
        xcen = (xmax+xmin)/2.0
        ycen = (ymax+ymin)/2.0
        ra_cnt, dec_cnt = telcoord.detpix2radec(xcen, ycen, attrot, teldef)
        
        
        #vdetx = np.array([xmin, xmax, xmax, xmin])
        #vdety = np.array([ymin, ymin, ymax, ymax])

        ### FOV bottom 
        v1detx = np.linspace(xmin, xmax, ngrid)
        v1dety = np.linspace(ymin, ymin, ngrid)
        ### FOV right 
        v2detx = np.linspace(xmax, xmax, ngrid)
        v2dety = np.linspace(ymin, ymax, ngrid)
        ### FOV top
        v3detx = np.linspace(xmax, xmin, ngrid)
        v3dety = np.linspace(ymax, ymax, ngrid)
        ### FOV left
        v4detx = np.linspace(xmin, xmin, ngrid)
        v4dety = np.linspace(ymax, ymin, ngrid)

        vdetx = np.array([v1detx,v2detx,v3detx,v4detx]).flatten()
        vdety = np.array([v1dety,v2dety,v3dety,v4dety]).flatten()
        npos = len(vdetx)

        vra, vdec = telcoord.detpix2radec(vdetx, vdety, attrot, teldef)

        if coordsys.upper()[:3]=='GAL' :
            #skycoord_fk5 = SkyCoord(vra, vdec, frame='fk5', unit='deg')
            skycoord_icrs = SkyCoord(vra, vdec, frame='icrs', unit='deg')
            skycoord_gal = skycoord_icrs.galactic
            vxcrd = skycoord_gal.l.degree
            vycrd = skycoord_gal.b.degree

            
            ### modify l range from (0,360) to (-180,180)
            vxcrd -= (180.<vxcrd)*360.
            
            skycoord_icrs = SkyCoord(ra_cnt, dec_cnt, frame='icrs', unit='deg')
            skycoord_gal = skycoord_icrs.galactic
            xcrd_cnt = skycoord_gal.l.degree
            ycrd_cnt = skycoord_gal.b.degree
            
            ### modify l range from (0,360) to (-180,180)
            xcrd_cnt -= (180.<xcrd_cnt)*360.

        else :
            vxcrd = vra
            vycrd = vdec
            xcrd_min = 0.
            xcrd_max = 360.

            xcrd_cnt = ra_cnt
            ycrd_cnt = dec_cnt


        ### detector ID
        regfile.write('text %lf %lf # text={%d} color=%s\n'%(xcrd_cnt, ycrd_cnt, teldef.detid, colorname))

        for ipos in range(npos) :
            ipos1 = ipos
            ipos2 = (ipos+1)%npos
            xcrd1 = vxcrd[ipos1]
            xcrd2 = vxcrd[ipos2]
            ycrd1 = vycrd[ipos1]
            ycrd2 = vycrd[ipos2]
            if 180.<abs(xcrd1-xcrd2) :
                continue
                
            regfile.write('line(%lf,%lf,%lf,%lf) # color=%s\n'%(xcrd1, ycrd1, xcrd2, ycrd2, colorname))
            
        idx+=1

    regfile.close()

    return
