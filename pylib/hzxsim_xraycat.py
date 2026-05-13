#!/usr/bin/env python

import os, sys
import numpy as np

import astropy.io.fits as pyfits

from hzxsim_point import HZXSimPoint

from xraycatalog import CatmxmonSources, Cat2rxsSources


class HZXSimXraycat (HZXSimPoint) :
    def __init__(self) :
        self.xraysource = None ### default
        
        return
    
    def simEventGen_XraysrcList(self, xraysrclist, tstart, tstop, cos_theta_limit=0.9721, detmm_margin=33.) :

        for xraysrc in xraysrclist :
            objectName = xraysrc["objectName"].replace(" ", "_")
            ra  = xraysrc["ra"]
            dec = xraysrc["dec"]
            vdet_obj = self.calcObjectDetCoord(ra, dec)
            if vdet_obj[2] < cos_theta_limit :
                continue

            xsmodel = xraysrc["xsmodel"]
            xsparam = xraysrc["xsparam"]
            if isinstance(xsparam, str) :
                xsparam = eval(xsparam)

            lcfname = None ### default
            if "lcfile" in xraysrc.keys() :
                lcfname = xraysrc["lcfile"]
                
            print(objectName, ra, dec, vdet_obj, xsmodel, xsparam)

            if not (len(xsmodel)>0 and len(xsparam)>0) :
                print("Warning: %s = %s does not have a spectral model."%(mxname, objectName))
                continue

            if self.xraysource == None :
                self.initXraysource(objectName, ra, dec, xsmodel=xsmodel, xsparam=xsparam)
            else :
                self.setObjectNameRADEC(objectName, ra, dec)
                self.setXspecModel(xsmodel, xsparam)

            if lcfname is not None :
                if os.path.isfile(lcfname) :
                    print("Set lightcurve data from %s"%(lcfname))
                    ### read LC file with defalut parameters 
                    self.readLcfile(lcfname, hdunum=1, timecol="TIME", ratecol="RATE", convfact=1.0, elo=0.1, ehi=10.0)
                else :
                    print("Warinig: lightcurve file %s does not exist !"%(lcfname))


            #if self.hzxpsf.version_number == 0 :
            #    self.simEventGen_psfv0(tstart, tstop, cos_theta_limit, detmm_margin)
            #else :  ### default
            self.simEventGen_bright(tstart, tstop, cos_theta_limit, detmm_margin)
        


    def simEventGen_MaximoniSources(self, mxmonSources, tstart, tstop, cos_theta_limit=0.9721, detmm_margin=33., mxlcdir="") :

        ####  MAXI catalog source  ####
        vvdet_obj = self.calcObjectDetCoord(mxmonSources.vra, mxmonSources.vdec)
        #print(vvdet_obj)
        
        num_obj = mxmonSources.nobj ## top num sources
        for idx in range(num_obj) :
            ### check target off-axis angle
            if vvdet_obj[2][idx] < cos_theta_limit :
                continue

            ra  = mxmonSources.vra[idx]
            dec = mxmonSources.vdec[idx]
            mxname = mxmonSources.vname[idx]
            objectName = mxmonSources.vidname[idx].replace(" ", "_")
            xsmodstr = mxmonSources.vxsmodel[idx]
            xsmodpar = mxmonSources.vxsmodpar[idx]

            if not (len(xsmodstr)>0 and len(xsmodpar)>0) :
                print("Warning: %s = %s does not have a spectral model."%(mxname, objectName))
                continue

            #print(xsmodstr, xsmodpar)
            if self.xraysource == None :
                self.initXraysource(objectName, ra, dec, xsmodel=xsmodstr, xsparam=eval(xsmodpar))
            else :
                self.setObjectNameRADEC(objectName, ra, dec)
                self.setXspecModel(xsmodstr, eval(xsmodpar))

                
            ### check MAXI lightcurve data
            lcbinsize = mxmonSources.vxslcbin[idx]
            lcdatfname = ""  ## default
            if mxlcdir!="" :
                if lcbinsize == 0.0 :
                    lcdatfname = mxlcdir+os.sep+"%s/mxglc_scan.fits"%(mxname)
                elif lcbinsize > 0.0 :
                    lcdatfname = mxlcdir+os.sep+"%s/mxglc_bin%.1lfh.fits"%(mxname, lcbinsize)

            #print("lcdatfname = ", lcdatfname)
            if len(lcdatfname)>0 :
                if os.path.isfile(lcdatfname) :
                    #elo = 2.0
                    #ehi = 4.0
                    #convfact = -1.0  ## calculate PhotonFlux-to-Rate factor in read LC function 
                    print("Set lightcurve data from %s"%(lcdatfname))
                    self.readLcfile(lcdatfname, hdunum=1, timecol='MJD', ratecol='RATE', convfact=-1.0, elo=2.0, ehi=4.0)

            #if self.hzxpsf.version_number == 0 :
            #    self.simEventGen_psfv0(tstart, tstop, cos_theta_limit, detmm_margin)
            #else :  ### default
            self.simEventGen_bright(tstart, tstop, cos_theta_limit, detmm_margin)
        
        
    def simEventGen_RosatSources(self, cat2rxsSources, tstart, tstop, cos_theta_limit=0.9721, detmm_margin1=33., detmm_margin2=3., mxmonSources=None):
        
            vvdet_obj = self.calcObjectDetCoord(cat2rxsSources.vra, cat2rxsSources.vdec)
            #print(vvdet_obj)

            num_obj = cat2rxsSources.nobj ## top num sources
            
            ### loop for catalog sources
            for idx in range(num_obj) :
                if vvdet_obj[2][idx] < cos_theta_limit :
                    continue
                
                ra  = cat2rxsSources.vra[idx]
                dec = cat2rxsSources.vdec[idx]
                rxsjname = cat2rxsSources.vname[idx].split()[1]
                objectName = cat2rxsSources.vname[idx].replace(" ","_")
                rate_hzx = cat2rxsSources.vrate_hzx[idx]

                if rate_hzx < 1.0e-3 : ### faint limit
                    continue
                
                ### check MAXI catalog
                if mxmonSources!=None :
                    mx_irow = mxmonSources.find_rxsobj(rxsjname, 2) # 2RXS
                    if mx_irow>=0 :
                        ### check spectral model with HZX rate
                        if  mxmonSources.vrate_hzx[mx_irow] > 0 :
                            print("2RXS %s exists in MAXI monitor catalog %s"%(rxsjname, mxmonSources.vname[mx_irow]))
                            continue
                

                if self.xraysource == None :
                    self.initXraysource(objectName, ra, dec) #, xsmodel=xsmodstr, xsparam=eval(xsmodpar))
                else :
                    self.setObjectNameRADEC(objectName, ra, dec)
                    #self.setXspecModel(xsmodstr, eval(xsmodpar))

                ### check if object is in FOV 
                if 0.02<rate_hzx : ### approximately top 1000 bright objects.
                    #this_detmm_margin=33.
                    detmm_margin = detmm_margin1
                else :
                    detmm_margin = detmm_margin2

                if not self.checkObjectFOV(cos_theta_limit, detmm_margin) :
                    print("Target object %s is out of FOV"%(objectName))
                    continue
                
                nh    = cat2rxsSources.vnh_pow[idx]
                gamma = cat2rxsSources.vgam_pow[idx]
                norm  = cat2rxsSources.vnorm_pow[idx]
                
                xsmodel = "wabs*power"
                xsparam = (nh, gamma, norm)
                self.setXspecModel(xsmodel, xsparam)

                #self.simEventGen_bright(tstart, tstop, cos_theta_limit, detmm_margin)
                #if self.hzxpsf.version_number == 0 :
                #    self.simEventGen_psfv0(tstart, tstop, cos_theta_limit, detmm_margin)
                #else :  ### default
                self.simEventGen_bright(tstart, tstop, cos_theta_limit, detmm_margin)

                
    
if __name__=="__main__" :


    #from hzxpsf import HZXPSF
    ###from hzxpsf_v0 import HZXPSF_v0 as HZXPSF
    #from hzxcaldb import Hzxcaldb
    #from hzxevtproc import HZXEvtProc
    
    from hzxpsf import HZXPSF, get_psffilelist
    from hzxcaldb import HzxCaldb
    from evtproc import EventProc

    hzxsimdir = os.getenv("HZXSIMDIR")

    ### MAXI monitor source catalog
    catfname = hzxsimdir+os.sep+"simdb/catmxmon.fits"
    catmxmonSources = CatmxmonSources(catfname)
    #mxlcdir_local = "/mnt/d1/mutsumi/maxi/catmxmon/v7.7l"
    #mxlcdir_local = "/Users/sugizaki/work/hiz/catmxmon/v7.7l"
    mxlcdir = hzxsimdir+os.sep+"simdb/mxgsclc"
    

    ### ROSAT 2RXS catalog
    catfname = hzxsimdir+os.sep+"simdb/cat2rxs_plparcut.fits"
    cat2rxsSources = Cat2rxsSources(catfname)


    #hzxsim = HZXSim() ### default
    hzxsim = HZXSimXraycat() ### default
    evtproc = EventProc() ### event processor

    ### init PSF database
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psffits"
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    psffitsdir =  hzxsimdir+os.sep+"refdata/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    psf_sigma = 0.2 # mmm
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)
    

    ###  Observation Time and Attitude
    #qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173] ### Crab region, detid=7 includes Crab  
    #qparam = [-0.82470631,  0.00189276,  0.00301748, -0.56555001]  #### Galactic center region
    qparam = [0.393458311142714, 0.139219718575607, 0.303128546265385, 0.856692191975324]  #### Cygnus region, detid=1 includes Cyg X-1
    tstart  = 0.0
    #tstop   = 1000.0
    tstop   = 10000.0
    objectName = "XRAYSKY"

    hzxsim.setAttQuat(qparam, inv=False)

    
    ### HZX simulate
    #detid_list = range(1,17)
    detid_list = [1]
    for detid in detid_list :

        detnam = "XM%02d"%(detid)
        print(detnam)

        ### init CALDB
        teldefname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
        rmffname   = hzxsimdir+os.sep+"refdata/hzx_erosita.rmf"
        arffname   = hzxsimdir+os.sep+ "refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
        vigfname   = hzxsimdir+os.sep+"refdata/hzx_vig_withframe_imgall.fits"
        hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
        #hzxcaldb  = HzxCaldb(detid)
        hzxsim.setCaldb(hzxcaldb)

        teldef = hzxcaldb.teldef
        
        ##hzxcaldb = Hzxcaldb(cmosid)
        ##hzxsim.setCaldb(hzxcaldb)
        #teldef = hzxcaldb.teldef
        backscale = hzxcaldb.teldef.det_xsiz*teldef.det_xscl * teldef.det_ysiz*teldef.det_yscl * 0.01 # (mm^2->cm^2) 

        
        ### init Event data
        hzxsim.initEventData()


        ### MAXI monitor source
        hzxsim.simEventGen_MaximoniSources(catmxmonSources, tstart, tstop, mxlcdir=mxlcdir)

        ### ROSAT source
        hzxsim.simEventGen_RosatSources(cat2rxsSources, tstart, tstop, mxmonSources=catmxmonSources)
        
        ### Detector response simulation
        hzxsim.simDetectorResp()

        ### Init data processing
        #if hzxsim.hzxpsf.version_number == 0 :
        #    outfname = "temp_psfv0_epw%d.evt"%(cmosid)
        #else :
        outfname = "xraycat_xm%02d.evt"%(detid)
        evtproc.setObsParams(teldef, tstart, tstop, qparam, objectName, inv=False)
        ### Store date to file
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)

