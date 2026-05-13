#!/usr/bin/env python

import os
import numpy as np
from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
import astropy.wcs as wcs

import xspec

#import datetime
#from astropy.time import Time
#import eptime
from mstime import MSTime, mstime_start, mstime_stop


from teldef import Teldef
import telcoord

from xraysource import XraySpectrum, LightCurve, XraySource


from hzxsim_base import HZXSim_base


class HZXSimPoint (HZXSim_base) :
    #def __init__(self, hzxcaldb=None, hzxpsfdb=None) :
    def __init__(self) :
        self.xraysource = None ### default
        return

    
    def initXraysource(self, objectName="TestObject", ra=0.0, dec=0.0, tstart=mstime_start, tstop=mstime_stop, xsmodel="phabs*power", xsparam=(0.3, 2.1, 1.0)) :
        if self.caldb == None :
            print("Error in HZXSimPount initXraysource. Caldb is not initizalized")
            return

        self.xraysource = XraySource(objectName, ra, dec, tstart, tstop, self.caldb.rmf, self.caldb.arf, xsmodel, xsparam)
        if self.attrot!=None :
            self.setObjectNameRADEC(objectName, ra, dec)

        
    def calcObjectDetCoord(self, vra, vdec) :
        vvdet_obj = telcoord.radec2vdet(vra, vdec, self.attrot, self.caldb.teldef.alignrot)
        return vvdet_obj

    
    def setObjectNameRADEC(self, objectName, ra, dec) :
        self.xraysource.setObjectNameRADEC(objectName, ra, dec)

        ### Calculate vdet_obj, detx_obj, dety_obj
        self.vdet_obj = telcoord.radec2vdet(ra, dec, self.attrot, self.caldb.teldef.alignrot)
        print(self.vdet_obj)
        
        self.detx_obj, self.dety_obj = telcoord.vdet2detmm(self.vdet_obj, self.caldb.teldef)
        print("target detxy (mm)  = %lf, %lf"%(self.detx_obj, self.dety_obj))
        self.detxpix_obj, self.detypix_obj = telcoord.detmm2detpix(self.detx_obj, self.dety_obj, self.caldb.teldef)
        print("target detxy (pixel)  = %lf, %lf"%(self.detxpix_obj, self.detypix_obj))

        return


    def setXspecModel(self, xsmodel, xsparam) :
        self.xraysource.setXspecModel(xsmodel, xsparam)
        return
    
    def readLcfile(self, xlcdatfname, hdunum=1, timecol="TIME", ratecol="RATE", convfact=1.0, elo=0.1, ehi=10.0) :
        self.xraysource.readLcfile(xlcdatfname, hdunum, timecol, ratecol, convfact, elo, ehi)
        return


    def checkObjectFOV(self, cos_theta_limit=0.9721, detmm_margin=33) :
        isObjectFOV=False ### default
        if self.vdet_obj[2] > cos_theta_limit : # vdet[0]=cos(theta)
            if (self.detxmi - detmm_margin) < self.detx_obj and self.detx_obj < (self.detxma + detmm_margin) and (self.detymi - detmm_margin) < self.dety_obj and self.dety_obj < (self.detyma + detmm_margin) :
                #print("Object: %s is out of PSF area"%(self.xraysource.objectName))
                isObjectFOV=True 
        return isObjectFOV
        
        
    
    def simEventGen_bright(self, tstart, tstop, cos_theta_limit=0.9721, detmm_margin=33.) : 
        ### cos_thete_limit = np.cos(np.arctan(64*np.sqrt(2)/375.) in default
        ### detmm_margin = 33 mm in default

        if not self.checkObjectFOV(cos_theta_limit, detmm_margin) :
            print("Target object is out of FOV")
            return
        
        ### get xraysource RA, Dec
        ra  = self.xraysource.ra
        dec = self.xraysource.dec

        #### initialize PSF 
        self.hzxpsf.initPSFevts(self.detx_obj, self.dety_obj)

        ### create empty array
        ### PSF event generateor
        for iband in range(self.hzxpsf.nbands) :
            elo = self.hzxpsf.elo_list[iband]
            ehi = self.hzxpsf.ehi_list[iband]
            psfareafact = self.hzxpsf.areafact_list[iband]
            if 0.0<psfareafact :
                this_vtime, this_vene = self.xraysource.genEvents(tstart, tstop, elo, ehi, psfareafact)
                nevts = len(this_vtime)

                if nevts > 0 :
                    ### PSF response
                    print("PSF response for nevts = %d"%(nevts))
                    this_vdetx, this_vdety = self.hzxpsf.genPSFevts_iband(iband, nevts)

                    self.vtime = np.append(self.vtime, this_vtime)
                    self.vene  = np.append(self.vene,  this_vene)
                    self.vdetx = np.append(self.vdetx, this_vdetx)
                    self.vdety = np.append(self.vdety, this_vdety)
                    self.vra   = np.append(self.vra,  np.zeros(nevts)+ra)
                    self.vdec  = np.append(self.vdec, np.zeros(nevts)+dec)
                    ### append empty vpha array
                    self.vpha  = np.append(self.vpha, np.zeros(nevts, 'int16'))

        return

    
    def simEventGen_faint(self, tstart, tstop, cos_theta_limit=0.99, detmm_margin=0.2, init_psf=True) : #, init_data=False) :
        ### cos_thete_limit = np.cos(np.arctan((32+5)*np.sqrt(2)/375.)

        ### check source incident angle
        if self.vdet_obj[2] < cos_theta_limit : # vdet[2]=cos(theta)
            print("Target object is out of FOV")
            return

        ### check source position in CMOS region
        if not ( (self.detxmi - detmm_margin) < self.detx_obj and self.detx_obj < (self.detxma + detmm_margin)
                 and (self.detymi - detmm_margin) < self.dety_obj and self.dety_obj < (self.detyma + detmm_margin) ) :
            print("Object: %s is out of imager area"%(self.xraysource.objectName))
            return
        
        ### Photon generator PSF paramters
        #elo = 0.1
        #ehi = 10.0
        elo = 0.1   ### RMF energy range
        ehi = 12.0

        #### initialize PSF 
        if init_psf :
            self.hzxpsf.initPSFevts(0.0, 0.0)

        ### initialize vignetting function
        #vigfunc = self.caldb.hzxvig.get_vigfunc(self.xraysource.spectrum.vene, self.detxpix_obj, self.detypix_obj)
        vigfunc = self.caldb.hzxvig.get_vigfunc(self.xraysource.spectrum.vene, self.detx_obj, self.dety_obj)
        print(vigfunc) ### test
        
        ### get xraysource RA, Dec
        ra  = self.xraysource.ra
        dec = self.xraysource.dec

        #temp_vtime, temp_vene = self.xraysource.genEvents(tstart, tstop, elo, ehi, psfareafact)
        temp_vtime, temp_vene = self.xraysource.genEvents(tstart, tstop, elo, ehi, vigfunc)
        temp_nevts = len(temp_vtime)
        if temp_nevts > 0 :

            temp_detx, temp_dety, temp_ngoodevts = self.hzxpsf.genPSFevts_defframe(temp_vene, self.detx_obj, self.dety_obj, self.detxmi, self.detxma, self.detymi, self.detyma)
            print("sim faint. nevts=%d  ngoodevts=%d"%(temp_nevts, temp_ngoodevts)) 
            if temp_ngoodevts>0 :
                self.vtime = np.append(self.vtime, temp_vtime)
                self.vene  = np.append(self.vene,  temp_vene)
                self.vdetx = np.append(self.vdetx, temp_detx)
                self.vdety = np.append(self.vdety, temp_dety)
                self.vra   = np.append(self.vra,  np.zeros(temp_nevts)+ra)
                self.vdec  = np.append(self.vdec, np.zeros(temp_nevts)+dec)
                ### append empty vpha array
                self.vpha  = np.append(self.vpha, np.zeros(temp_nevts, 'int16'))
        return


    

if __name__=="__main__" :


    from hzxpsf import HZXPSF, get_psffilelist
    from hzxcaldb import HzxCaldb
    from evtproc import EventProc

    ### default simulator
    hzxsim = HZXSimPoint() 

    ### Event process
    evtproc = EventProc()


    ###  Observation Time and Attitude
    detid  = 7
    detnam = "XM%02d"%(detid)
    qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173]

    tstart  = 0.0
    tstop   = 1000.0
    
    ### Observation target
    objectName = "Crab"
    ra  = 83.633
    dec = 22.014
    

    hzxsimdir = os.getenv("HZXSIMDIR")
    ### init PSF database
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psffits"
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    psffitsdir = hzxsimdir+os.sep+"refdata/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    psf_sigma = 0.2 # mmm
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)

    ### init CALDB
    #teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
    #rmffname   = "../refdata/hzx_erosita.rmf"
    #arffname   = "../refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
    #vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"

    teldefname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
    rmffname   = hzxsimdir+os.sep+"refdata/hzx_erosita.rmf"
    arffname   = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
    vigfname   = hzxsimdir+os.sep+"refdata/hzx_vig_withframe_imgall.fits"

    hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
    #hzxcaldb  = HzxCaldb(detid)
    hzxsim.setCaldb(hzxcaldb)


    ### inv=True for attiude data
    ### set inv=False for the old qparam
    #hzxsim.setAttQuat(qparam, inv=False)  
    hzxsim.setAttQuat(qparam, inv=False)  

    
    ### set target position, X-ray spectrum, light curve
    #hzxsim.setObjectNameRADEC(objectName, ra, dec)
    #hzxsim.setXspecModel("phabs*power", (0.3, 2.1, 10.0))
    hzxsim.initXraysource(objectName, ra, dec, xsmodel="phabs*power", xsparam=(0.3, 2.1, 10.0))



    ### Faint source case
    outfname = "crab_faint.evt"
    if os.path.isfile(outfname) :
        print("Output file=%s already exists"%(outfname))
    else :
        hzxsim.initEventData()
        hzxsim.simEventGen_faint(tstart, tstop)
        hzxsim.simDetectorResp()
        ### process
        evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, inv=False)
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)


    ### Bright source case
    outfname = "crab_bright.evt"
    if os.path.isfile(outfname) :
        print("Output file=%s already exists"%(outfname))
    else :
        hzxsim.initEventData()
        hzxsim.simEventGen_bright(tstart, tstop)
        hzxsim.simDetectorResp()
        ## process
        evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, inv=False)
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)
    

    
