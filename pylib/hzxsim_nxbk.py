#!/usr/bin/env python

import os
import numpy as np
#from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
#import astropy.wcs as wcs


import telcoord
from hzxsim_base import HZXSim_base


class HZXSimNXBK (HZXSim_base) :
    def __init__(self) :
        self.vratespec = np.array([])
        self.vratecdf  = np.array([])
        return
    
    
    def readNXBKSpec(self, specfname, backscale=1.0) :

        ### open NXB spectrum file
        hdul = pyfits.open(specfname)
        is_ebounds = False   ### default
        for hdu in hdul :
            if "EXTNAME" in hdu.header.keys():
                if hdu.header["EXTNAME"]=="EBOUNDS" :
                    is_ebounds = True
                
        if is_ebounds==False :  ### default
            print("NXB file is a default spectrum.")
            hdu = hdul["SPECTRUM"] ### Spectrum HDU
            nbins = hdu.header["NAXIS2"]
            exposure = hdu.header["EXPOSURE"]
            vchan = hdu.data.field("CHANNEL").astype('i')
            vrate = hdu.data.field("COUNTS").astype('d')/exposure

            if np.all(self.caldb.rmf.vchan == vchan) : ### Chennel must be exactly same.
                self.vratespec = vrate * backscale
                self.vratecdf = np.cumsum(self.vratespec)
                self.vratecdf/=self.vratecdf[-1]
            else :
                print("Error in readNXBKSpec. inconsistent channnels between specfile and rmf")
                self.vratespec = np.array([])
                self.vratecdf  = np.array([])

        else :
            print("NXB file has a EBOUNDS extension")
            hdu = hdul["EBOUNDS"]
            nchans = hdu.header["NAXIS2"]
            vemin = hdu.data.field("E_MIN").astype('d')
            vemax = hdu.data.field("E_MAX").astype('d')
            vecen = (vemin+vemax)/2.0

            hdu = hdul["SPECTRUM"]
            exposure = hdu.header["EXPOSURE"]
            vrate = hdu.data.field("COUNTS").astype('d')/exposure

            self.vratespec = np.zeros(self.caldb.rmf.nchans)
            for ichan in range(nchans) :
                ecen = vecen[ichan]
                ichan_resp = np.searchsorted(self.caldb.rmf.vchan_emax, ecen)
                #print(ecen, ichan_resp)
                self.vratespec[ichan_resp] += vrate[ichan]

            self.vratespec *= backscale
            self.vratecdf = np.cumsum(self.vratespec)
            self.vratecdf/=self.vratecdf[-1]
            
            
        hdul.close()

        return
        
    
    def simEventGen_NXBK(self, tstart, tstop) :

        exposure = tstop - tstart

        nevts_expect = self.vratespec.sum() * exposure #* detarea
        nevts = np.random.poisson(nevts_expect)
        print("nevts_expect=%lf, nevts=%d"%(nevts_expect, nevts))

        #vrand = np.random.rand(nevts) ### fix bug in the correlation among vtime, vpha, vdetxpix
        vtime = tstart + np.random.rand(nevts) * exposure 
        vpha  = np.searchsorted(self.vratecdf, np.random.rand(nevts)).astype('int16')
        vdetxpix = (self.caldb.teldef.detxpix1 + np.random.rand(nevts) * self.caldb.teldef.det_xsiz).astype('i')
        #vrand2   = np.random.rand(nevts)  ### random numbers for dety 
        vdetypix = (self.caldb.teldef.detypix1 + np.random.rand(nevts) * self.caldb.teldef.det_ysiz).astype('i')
        vdetx, vdety = telcoord.detpix2detmm(vdetxpix, vdetypix, self.caldb.teldef)

        self.vtime = np.append(self.vtime, vtime)
        self.vpha  = np.append(self.vpha, vpha)
        #self.vdetxpix = np.append(self.vdetxpix, vdetxpix)
        #self.vdetypix = np.append(self.vdetypix, vdetypix)
        
        self.vene  = np.append(self.vene, np.zeros(nevts)-1.0)
        self.vra   = np.append(self.vra,  np.zeros(nevts)-999.0)
        self.vdec  = np.append(self.vdec, np.zeros(nevts)-999.0)
        self.vdetx = np.append(self.vdetx, vdetx)
        self.vdety = np.append(self.vdety, vdety)

        

        
if __name__=="__main__" :

    #from hzxpsf import HZXPSF, get_psffilelist
    from hzxcaldb import HzxCaldb
    from evtproc import EventProc


    hzxsim = HZXSimNXBK() ### default
    evtproc = EventProc() ### event processor
    
    ### HZX simulate
    objectName = "NXB"

    ###  Observation Time and Attitude
    tstart  = 0.0
    tstop   = 1000.0
    qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173]

    nxbspecfname = "../refdata/fake_bkg_exp1e6s.pha"
    nxbscale = 1.0

    ### init PSF database
    ##psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psffits"
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    #ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    #this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    #psf_sigma = 0.2 # mmm
    #hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    #hzxsim.setHzxpsf(hzxpsfdb)

    hzxsim.setAttQuat(qparam)  ### dummy
    
    detid_list  = [1,2]

    for detid in detid_list :
        detnam = "XM%02d"%(detid)
        print(detnam)

        ### init CALDB
        teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
        rmffname   = "../refdata/hzx_erosita.rmf"
        arffname   = "../refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
        vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"
        hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
        #hzxcaldb  = HzxCaldb(detid)
        hzxsim.setCaldb(hzxcaldb)


        hzxsim.initEventData()

        hzxsim.readNXBKSpec(nxbspecfname, nxbscale)
        hzxsim.simEventGen_NXBK(tstart, tstop) 


        #if len(vtime)>0 :
        outfname = "temp_nxbkg_xm%02d.evt"%(detid)
        if os.path.isfile(outfname) :
            os.remove(outfname)


        ### Detector response simulation
        hzxsim.simDetectorResp()

        ### process
        evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, inv=False)
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)



