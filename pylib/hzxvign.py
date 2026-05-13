#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits


class HZXVign :
    def __init__(self, vigfname, eband_le=0.1, eband_he=12.0) :

        self.eband_le = eband_le
        self.eband_he = eband_he

        #self.nene = len(elist)
        #self.vene = np.array(elist)
        ### check vignetting image format

        hdul = pyfits.open(vigfname)

        ### get energy list
        hdu = hdul[1]
        self.nene = hdu.header["NAXIS2"]
        self.vene = hdu.data.field("Energy")

        ### get image arr
        hdu = hdul[0]
        self.imgarr = hdu.data
        
        ### need to check X,Y order
        self.ny = hdu.header["NAXIS1"]  ### axis1 = Y
        self.nx = hdu.header["NAXIS2"]  ### axis2 = X
        self.crpix1 = hdu.header["CRPIX1"] 
        self.crpix2 = hdu.header["CRPIX2"] 
        self.crval1 = hdu.header["CRVAL1"] 
        self.crval2 = hdu.header["CRVAL2"] 
        self.cdelt1 = hdu.header["CDELT1"] 
        self.cdelt2 = hdu.header["CDELT2"] 

        hdul.close()

        ### normalize vignetting data
        for iene in range(self.nene) :
            self.imgarr[iene]/=self.imgarr[iene].max()
            
        

        self.elo_list = np.zeros(self.nene)
        self.ehi_list = np.zeros(self.nene)
        
        for iband in range(self.nene) :
            if iband==0 :
                self.elo_list[iband] = self.eband_le
            else :
                self.elo_list[iband] = (self.vene[iband-1] + self.vene[iband])/2.0
                
            if iband == self.nene -1 :
                self.ehi_list[iband] = self.eband_he
            else :
                self.ehi_list[iband] = (self.vene[iband] + self.vene[iband+1])/2.0

        return

    
    def get_vigfact(self, ene, xpix, ypix) :
        ### xpix = detx in mm, ypix = DETY in mm
        
        ### interpolate energy bin
        iz1 = np.argmin(np.fabs(self.vene-ene))
        ene1 = self.vene[iz1]
        if ene==ene1 :
            iz2 = -1
        elif ene < ene1 :
            iz2 = iz1-1
        else :
            iz2 = iz1+1

        if iz2<0 or self.nene-1<=iz2 :
            iz2 = -1
            ene2 = -1.
        else :
            ene2 = self.vene[iz2]

        #if not (ene1 <= ene and ene <= ene2) :
        #print(ene1, ene2, ene)

        ### interpolate image bin
        iy1 = int(np.floor( (ypix - self.crval1)/self.cdelt1 + self.crpix1 ))
        iy2 = iy1 + 1

        ix1 = int(np.floor( (xpix - self.crval2)/self.cdelt2 + self.crpix2 ))
        ix2 = ix1 + 1

        if 0<=iy1 and iy2<self.ny and 0<=ix1 and ix2<self.nx :

            ypix1 = (iy1 - self.crpix1)*self.cdelt1 + self.crval1
            ypix2 = (iy2 - self.crpix1)*self.cdelt1 + self.crval1
            
            xpix1 = (ix1 - self.crpix2)*self.cdelt2 + self.crval2
            xpix2 = (ix2 - self.crpix2)*self.cdelt2 + self.crval2
            #if not (xpix1 <= xpix and xpix <= xpix2) :
            #    print(ix1, ix2, xpix1, xpix2, xpix)
            
            z1yx11 = self.imgarr[iz1][iy1][ix1]
            z1yx12 = self.imgarr[iz1][iy1][ix2]
            z1y1 = (z1yx12-z1yx11)/(xpix2-xpix1)*(xpix-xpix1) + z1yx11

            z1yx21 = self.imgarr[iz1][iy2][ix1]
            z1yx22 = self.imgarr[iz1][iy2][ix2]
            z1y2 = (z1yx22-z1yx21)/(xpix2-xpix1)*(xpix-xpix1) + z1yx21

            z1 = (z1y2-z1y1)/(ypix2-ypix1)*(ypix-ypix1) + z1y1

            if iz2 < 0 :
                zval = z1

            else :
                z2yx11 = self.imgarr[iz2][iy1][ix1]
                z2yx12 = self.imgarr[iz2][iy1][ix2]
                z2y1 = (z2yx12-z2yx11)/(xpix2-xpix1)*(xpix-xpix1) + z2yx11

                z2yx21 = self.imgarr[iz2][iy2][ix1]
                z2yx22 = self.imgarr[iz2][iy2][ix2]
                z2y2 = (z2yx22-z2yx21)/(xpix2-xpix1)*(xpix-xpix1) + z2yx21

                z2 = (z2y2-z2y1)/(ypix2-ypix1)*(ypix-ypix1) + z2y1
                zval = (z2-z1)/(ene2-ene1)*(ene-ene1) + z1

                #if zval<0.0 :
                #    #print(z1y1, z1y2, z2y1, z2y2)
                #    print(z1y1, z1yx11, z1yx12, xpix1, xpix2, xpix)

        else :
            zval = 0.0

        return zval

    
    def get_vigfunc(self, vene, xpix, ypix) :
        ## xpix, ypin in mm
        vfact = np.zeros(self.nene)
        for iene in range(self.nene) :
            vfact[iene] = self.get_vigfact(self.vene[iene], xpix, ypix)

        vfunc = np.interp(vene, self.vene, vfact)
        return vfunc

    
if __name__=="__main__" :

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    #vigfname = "../refdata/hizx_vig_imgall.fits"
    vigfname = "../refdata/hzx_vig_withframe_imgall.fits"
    vign = HZXVign(vigfname)

    for idx in range(vign.nene) :
        print("%2d: %5.2lf  %5.2lf"%(idx, vign.elo_list[idx], vign.ehi_list[idx]))
    

    ### 
    nx = 500
    ny = 500
    xmi = -35
    xma =  35
    ymi = -35
    yma =  35
    
    vobjy = np.linspace(xmi, xma, nx)
    vobjx = np.linspace(ymi, yma, nx)

    img = np.zeros((ny, nx))

    #ekev = 0.5
    ekev = 1.0
    #ekev = 2.0
    for iy in range(ny) :
        objy = vobjy[iy]
        for ix in range(nx) :
            objx = vobjy[ix]
            vigfact = vign.get_vigfact(ekev, objx, objy) 
            img[iy][ix] = vigfact

    fig, ax = plt.subplots(figsize=(5,5), subplot_kw={"projection": "3d"})
    xx, yy = np.meshgrid(vobjx, vobjy)
    ax.plot_surface(xx, yy, img, cmap='viridis')
        
    plt.ion()
    plt.show()

