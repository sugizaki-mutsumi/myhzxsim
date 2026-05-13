#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits

from astropy.coordinates import SkyCoord


class CatmxmonSources :
    def __init__(self, catfname, sort_hzxrate=True) :

        fptr = pyfits.open(catfname)
        hdu = fptr[1]
        nrows = hdu.header["NAXIS2"] ### default

        ###  expected HZX observation rate and photon flux 
        vrate_hzx = hdu.data.field("RATE_OBS").astype("d")
        vpflux    = hdu.data.field("PFLUX_OBS").astype("d")

        ### descending order 
        if sort_hzxrate :
            idxarr = np.flip(np.argsort(vrate_hzx)) 
        else :
            idxarr = np.flip(np.argsort(vpflux)) 
            
        self.vrate_hzx = vrate_hzx[idxarr]
        self.vpflux0 = hdu.data.field("PFLUX0").astype("d")[idxarr]
        self.vpflux1 = hdu.data.field("PFLUX1").astype("d")[idxarr]
        self.vpflux2 = hdu.data.field("PFLUX2").astype("d")[idxarr]
        self.vpflux3 = hdu.data.field("PFLUX3").astype("d")[idxarr]


        self.vname = hdu.data.field("MXNAME")[idxarr]
        self.vra     = hdu.data.field("RA").astype("d")[idxarr]
        self.vdec    = hdu.data.field("DEC").astype("d")[idxarr]

        self.vjname_1rxs = hdu.data.field("NAME_1RXS")[idxarr]
        self.vrate_1rxs  = hdu.data.field("RATE_1RXS").astype("d")[idxarr]

        self.vjname_2rxs = hdu.data.field("NAME_2RXS")[idxarr]
        self.vrate_2rxs  = hdu.data.field("RATE_2RXS").astype("d")[idxarr]

        self.vxsmodel  = hdu.data.field("XSMODEL")[idxarr]
        self.vxsmodpar = hdu.data.field("XSMODPAR")[idxarr]
        self.vxslcbin  = hdu.data.field("XSLCBIN").astype("d")[idxarr]

        self.vidname  = hdu.data.field("IDNAME")[idxarr]

        fptr.close()

        self.nobj = nrows
        self.vcoord = SkyCoord(self.vra, self.vdec, frame='fk5', unit='deg')
        return


    def find_jname_id(self, jname) :
        idxlist = np.where(self.vname==jname) 
        if len(idxlist[0])==0 :
            irow = -1
        else :
            irow = idxlist[0][0]
        return irow

        
    def find_rxsobj(self, rxsjname, rxsid=1) :
        if rxsid==1 :
            vjname_rxs = self.vjname_1rxs
        else :
            vjname_rxs = self.vjname_2rxs

        idxlist = np.where(vjname_rxs==rxsjname) 

        if len(idxlist[0])==0 :
            irow = -1
        else :
            irow = idxlist[0][0]
        return irow
            
    
    def get_closest_id(self, ra, dec) :
        this_coord = SkyCoord(ra, dec, frame='fk5', unit='deg')
        vsep = this_coord.separation(self.vcoord)
        irow = vsep.argmin()
        sepdeg  = vsep[irow].deg
        return irow, sepdeg


    
        

class Cat2rxsSources :
    def __init__(self, catfname, sort_hzxrate=True) :

        fptr = pyfits.open(catfname)
        hdu = fptr[1]
        nrows = hdu.header["NAXIS2"] ### deafault

        vrate = hdu.data.field("RATE").astype("d")
        ### rate obs
        vrate_hzx = hdu.data.field("RATE1").astype("d")  ### HZX Obs rate

        ### descending order 
        if sort_hzxrate :
            idxarr = np.flip(np.argsort(vrate_hzx)) #[idx_start:idx_stop] ### desending order
        else :
            idxarr = np.flip(np.argsort(vrate)) #[idx_start:idx_stop] ### desending order
            
        self.vrate     = vrate[idxarr]
        self.vrate_hzx = vrate_hzx[idxarr]
        self.vname = hdu.data.field("IAU_NAME")[idxarr]
        self.vra   = hdu.data.field("RA_DEG").astype("d")[idxarr]
        self.vdec  = hdu.data.field("DEC_DEG").astype("d")[idxarr]

        self.vnh_pow    = hdu.data.field("NH_fit_p").astype("d")[idxarr]
        self.vnhe_pow   = hdu.data.field("NH_err_p").astype("d")[idxarr]
        self.vgam_pow   = hdu.data.field("GAMMA_p").astype("d")[idxarr]
        self.vgame_pow  = hdu.data.field("GAMMA_err_p").astype("d")[idxarr]
        self.vnorm_pow  = hdu.data.field("NORM_p").astype("d")[idxarr]
        self.vnorme_pow = hdu.data.field("NORM_err_p").astype("d")[idxarr]

        fptr.close()

        self.nobj = nrows
        self.vcoord = SkyCoord(self.vra, self.vdec, frame='fk5', unit='deg')
        return


    def find_jname_id(self, jname) :
        rxname = "2RXS %s"%(jname)
        idxlist = np.where(self.vname==rxname) 
        if len(idxlist[0])==0 :
            irow = -1
        else :
            irow = idxlist[0][0]
        return irow


    def get_closest_id(self, ra, dec) :
        this_coord = SkyCoord(ra, dec, frame='fk5', unit='deg')
        vsep = this_coord.separation(self.vcoord)
        irow = vsep.argmin()
        sepdeg  = vsep[irow].deg
        return irow, sepdeg


    
if __name__=="__main__" :


    #xtsimdir = os.getenv("HZXSIMDIR")
    simdir = ".."
    
    ### MAXI catalog
    catfname = simdir+os.sep+"simdb/catmxmon.fits"
    mxmonSources = CatmxmonSources(catfname)

    ### ROSAT 2RXS catalog
    catfname = simdir+os.sep+"simdb/cat2rxs_plparcut.fits"
    cat2rxsSources = Cat2rxsSources(catfname)


    num = 100
    ### MAXI
    print("MAXI monitor catalog, ra, dec, hzxrate_mx, hzxrate_2rxs, pflux0, rate_2rxs, rate_1rxs, idname")
    for idx in range(num) :
        srcjname = mxmonSources.vname[idx]
        idname = mxmonSources.vidname[idx]
        ra  = mxmonSources.vra[idx]
        dec = mxmonSources.vdec[idx]
        hzxrate_mx  = mxmonSources.vrate_hzx[idx]
        pflux0  = mxmonSources.vpflux0[idx]
        rate_1rxs = mxmonSources.vrate_1rxs[idx]
        rate_2rxs = mxmonSources.vrate_2rxs[idx]

        jname_2rxs = mxmonSources.vjname_2rxs[idx]
        idx_2rxs = cat2rxsSources.find_jname_id(jname_2rxs)
        if idx_2rxs>=0 :
            hzxrate_2rxs = cat2rxsSources.vrate_hzx[idx_2rxs]
        else :
            hzxrate_2rxs = -1.0

        print("%3d %9s %16s  %8.4lf %8.4lf  %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf  %s"%(idx, srcjname, jname_2rxs, ra, dec, hzxrate_mx, hzxrate_2rxs, pflux0, rate_2rxs, rate_1rxs, idname))
        
        
    #### ROSAT
    print("")
    print("ROSAT 2RXS catalog, ra, dec, hzxrate_mx, hzxrate_2rxs, rate_2rxs, idname")

    for idx in range(num) :
        srcjname = cat2rxsSources.vname[idx].split()[1]
        ra  = cat2rxsSources.vra[idx]
        dec = cat2rxsSources.vdec[idx]
        rate_2rxs    = cat2rxsSources.vrate[idx]
        hzxrate_2rxs = cat2rxsSources.vrate_hzx[idx]

        idx_mx = mxmonSources.find_rxsobj(srcjname, 2)
        if idx_mx>=0 :
            jname_mx   = mxmonSources.vname[idx_mx]
            hzxrate_mx = mxmonSources.vrate_hzx[idx_mx]
            idname = mxmonSources.vidname[idx_mx]
        else :
            jname_mx = "---"
            idname = "---"
            hzxrate_mx = -1.
        
            
        print("%3d %16s %9s  %8.4lf %8.4lf %8.3lf %8.3lf %8.3lf  %s"%(idx, srcjname, jname_mx, ra, dec, hzxrate_mx, hzxrate_2rxs, rate_2rxs, idname))
        

