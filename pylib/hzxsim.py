#!/usr/bin/env python

import os #, sys
import numpy as np
#from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits

from hzxsim_xraycat import HZXSimXraycat
from hzxsim_xrbk import HZXSimXRBK
from hzxsim_nxbk import HZXSimNXBK


class HZXSim (HZXSimXraycat, HZXSimXRBK, HZXSimNXBK) :
    def __init__(self) :
        self.xraysource = None ### default
        self.xrbkspec_list = []
        
        return
    
        
if __name__=="__main__" :

    from hzxpsf import HZXPSF, get_psffilelist
    from hzxcaldb import HzxCaldb
    from evtproc import EventProc
    
    hzxsim  = HZXSim() ### default
    evtproc = EventProc() ### event processor
    
    ### HZX simulate
    objectName = "XRAYSKY"

    ### init PSF database
    psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    psf_sigma = 0.2 # mmm
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)
    

    qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173] ### Crab region
    #qparam = [-0.82470631,  0.00189276,  0.00301748, -0.56555001]  #### Galactic center region
    tstart  = 0.0
    tstop   = 1000.0


    ### set Attitude
    #hzxsim.setAttQuat(qparam, inv=False) ### default: inv=False
    hzxsim.setAttQuat(qparam)
    
    
    ### set point source
    objectName = "Crab"
    ra  = 83.633
    dec = 22.014

    
    #hzxsimdir = os.getcwd()
    hzxsimdir = os.getenv("HZXSIMDIR")
    simdbdir = hzxsimdir+os.sep+"simdb"
    refdatadir = hzxsimdir+os.sep+"refdata"

    ### X-ray background emission data    
    #hpspecfname  = simdbdir+os.sep+"xrbkg_rayspec.fits"  
    hpspecfname  = simdbdir+os.sep+"xrbkg_rayspec_erobin.fits"  

    ### Non-X-ray background spectrum    
    nxbspecfname = refdatadir+os.sep+"fake_bkg_exp1e6s.pha"

    ### response file default
    rmffname   = refdatadir+os.sep+"hzx_erosita.rmf"
    ### arf file default
    arffname   = refdatadir+os.sep+"hzxarf_withframe_rn0.0nm_imgall.arf"
    ### vignetting data default
    vigfname   = refdatadir+os.sep+"hzx_vig_withframe_imgall.fits"

    teldefver = "20251216v002"

    #detid_list = range(1,2)
    detid_list = range(2,17)
    #cmosid_list = range(42,43)
    for detid in detid_list :

        print("XM%02d"%(detid))

        ### init CALDB
        teldefname = "../refdata/teldef/hiz_xm%02d_teldef_%s.fits"%(detid, teldefver)
        hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
        hzxsim.setCaldb(hzxcaldb)

        
        ### init Event data
        hzxsim.initEventData()


        ### point source events
        hzxsim.initXraysource(objectName, ra, dec, xsmodel="phabs*power", xsparam=(0.3, 2.1, 10.0))
        hzxsim.simEventGen_bright(tstart, tstop)
        #hzxsim.simDetectorResp()
        

        ### X-ray background events (from ROSAT all-sky survey data) 
        hzxsim.initXRBKSpec(hpspecfname)
        hzxsim.simEventGen_XRBK(tstart, tstop) 


        ### Non-X-ray background events (from simulation)
        hzxsim.readNXBKSpec(nxbspecfname, backscale=1)
        hzxsim.simEventGen_NXBK(tstart, tstop) 


        ### Detector response simulation
        hzxsim.simDetectorResp()

        ### Init data processing
        outfname = "temp_xm%02d.evt"%(detid)
        evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName) #, inv)
        ### Store date to file
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)

