#!/usr/bin/env python

import os
import numpy as np
from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
import astropy.wcs as wcs

import xspec

###import datetime
from astropy.time import Time
from mstime import MSTime, mstime_start, mstime_stop

from teldef import Teldef
import telcoord


from hzxpsf import HZXPSF, get_psffilelist
from hzxcaldb import HzxCaldb


class HZXSim_base :
    def __init__(self) :
    
        self.caldb = None ### default
        #if hzxcaldb!=None :
        #    ### initizalize caldb 
        #    self.setCaldb(hzxcaldb)
        
        self.hzxpsf = None ### default
        #if hzxpsfdb!=None :
        #    self.hzxpsf = hzxpsfdb
        

        ### initialize for each Source simulation ###
        self.attrot = None ### default


        ### initilize event data array
        self.initEventData()
        
        ##self.xraysource = None ### default

        return

    
    def setHzxpsf(self, hzxpsfdb) :
        self.hzxpsf = hzxpsfdb
        return
    

    def setCaldb(self, hzxcaldb) :

        self.caldb = hzxcaldb
        ### set imager DETXY range
        teldef = self.caldb.teldef 
        self.detxmi = teldef.det_xscl * (teldef.detxpix_min - teldef.detxpix_cen) 
        self.detxma = teldef.det_xscl * (teldef.detxpix_max - teldef.detxpix_cen) 
        self.detymi = teldef.det_yscl * (teldef.detypix_min - teldef.detypix_cen) 
        self.detyma = teldef.det_yscl * (teldef.detypix_max - teldef.detypix_cen) 

        #if self.xraysource!=None :
        #    self.xraysource.resetSpecresp(self.caldb.rmf, self.caldb.arf, rmfnorm=1)
        
        return

        
    #def setAttQuat(self, qparam, inv=True) :
    ### change default inv to False
    def setAttQuat(self, qparam, inv=False) :
        if inv :
            self.attrot = R.from_quat(qparam).inv()
        else :
            self.attrot = R.from_quat(qparam)
        
    def initEventData(self) :
        ### create empty array
        self.vra   = np.zeros(0, 'd')
        self.vdec  = np.zeros(0, 'd')
        self.vtime = np.zeros(0, 'd')
        self.vene  = np.zeros(0, 'd')
        self.vdetx = np.zeros(0, 'd')
        self.vdety = np.zeros(0, 'd')

        self.vpha     =  np.zeros(0,'i')
        self.vdetxpix =  np.zeros(0,'i')
        self.vdetypix =  np.zeros(0,'i')
        return        

    

    def simDetectorResp(self) : 

        ### check event data array
        # ...
        
        
        ### simulate imaging detector
        nevts = len(self.vtime)
        if nevts>0 :
            ### RMF response
            print("RMF response for nevts = %d"%(nevts))
            #self.vpha = self.rspmat.procEvents(self.vene)
            rmf = self.caldb.rmf

            ### modify later to use numpy function
            #self.vpha = np.zeros(nevts, 'int16')
            vrand =  np.random.rand(nevts)

            #viene = ((self.vene - rmf.elo)*rmf.ne / (rmf.ehi - rmf.elo)).astype('i')
            #np.clip(viene, 0, rmf.ne-1, None) ### clip lower/upper range of (0,1023)
            #for idx in range(nevts) :
            #    if 0.0<self.vene[idx] : ### exclude NXB data
            #        iene = viene[idx]
            #        ichan = np.searchsorted(rmf.respcdf[iene], vrand[idx])
            #        self.vpha[idx] = ichan
            
            #viene = ((self.vene - rmf.elo)*rmf.ne / (rmf.ehi - rmf.elo)).astype('i')
            #np.clip(viene, 0, rmf.ne-1, None) ### clip lower/upper range of (0,1023)
            viene = np.digitize(self.vene, rmf.veh) 
            for idx in range(nevts) :
                if 0.0<self.vene[idx] : ### exclude NXB data
                    iene = viene[idx]
                    ichan = np.searchsorted(rmf.respcdf[iene], vrand[idx])
                    self.vpha[idx] = ichan
            
            ### calculate vdetpix
            self.vdetxpix, self.vdetypix = telcoord.detmm2detpix(self.vdetx, self.vdety, self.caldb.teldef)

        else :
            self.vpha     =  np.zeros(0,'i')
            self.vdetxpix =  np.zeros(0,'i')
            self.vdetypix =  np.zeros(0,'i')
        return

    
    

if __name__=="__main__" :


    ### default simulator
    hzxsim = HZXSim_base() 

    ### Event process
    #evtproc = EventProc()


    detid  = 7
    detnam = "XM%02d"%(detid)
    qparam = [-0.70955226, -0.00353587,  0.00138493,  0.70464259]
    
    #tstart  = 103970520.0000197
    #tstop   = 103971425.9500197

    tstart  = 0.0
    tstop   = 1000.0

    
    ### Observation target
    objectName = "Crab"
    ra  = 83.633
    dec = 22.014
    

    ### init PSF database
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psffits"
    psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    psf_sigma = 0.2 # mmm
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)

    ### init CALDB
    teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
    rmffname   = "../refdata/hzx_erosita.rmf"
    arffname   = "../refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
    vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"
    hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
    #hzxcaldb  = HzxCaldb(detid)
    hzxsim.setCaldb(hzxcaldb)


    
    
