#!/usr/bin/env python

import os #, sys
import numpy as np
#from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
#import astropy.wcs as wcs

import healpy
#from astropy_healpix import healpy

#import xspec
#from astropy.time import Time
#import eptime

#from teldef import Teldef

import telcoord
from xraysource import XraySpectrum
from hzxsim_base import HZXSim_base


class HZXSimXRBK (HZXSim_base) :
    def __init__(self) :

        self.xrbkspec_list = []


    def init_xrbkspec_list(self) :
        ### check caldb
        if self.caldb == None :
            print("HZXSimXRBKG::init_xrbkspec_list error. Caldb has not been set")
            return 0

        ### check psf 
        if self.hzxpsf == None :
            print("HZXSimXRBKG::init_xrbkspec_list error. HZXPSF has not been set")
            return 0

        if self.hzxpsf.version_number == 0 :
            print("HZXSimXRBKG::init_xrbkspec_list error. HZXPSF version has to be >0")
            return 0

        
        teldef = self.caldb.teldef
        
        
        ####  setup BGD spatial parameters ####
        ####       common in all CMOS      ####
        self.nx = 24 #16
        self.ny = 24 #16
        #self.xmi = -0.015*2048  *1.5
        #self.xma =  0.015*2048  *1.5  
        #self.ymi = -0.015*2048  *1.5
        #self.yma =  0.015*2048  *1.5
        self.xmi = -teldef.det_xsiz * teldef.det_xscl *0.75
        self.xma =  teldef.det_xsiz * teldef.det_xscl *0.75
        self.ymi = -teldef.det_ysiz * teldef.det_yscl *0.75
        self.yma =  teldef.det_ysiz * teldef.det_yscl *0.75

        self.deltx = (self.xma-self.xmi)/self.nx
        self.delty = (self.yma-self.ymi)/self.ny

        self.pntx_list = np.zeros(self.nx * self.ny)
        self.pnty_list = np.zeros(self.nx * self.ny)

        #self.solangle = np.zeros(self.nx * self.ny)
        self.solidangle = np.degrees(self.deltx/teldef.focallen)*np.degrees(self.delty/teldef.focallen)*3600.  ### degre^2 -> arcmin^2

        self.nxy = self.nx * self.ny
        
        for iy in range(self.ny) :
            pnty = self.ymi + self.delty*(iy+0.5)   ### mm
            for ix in range(self.nx) :
                pntx = self.xmi + self.deltx*(ix+0.5)  ### mm
                ibin = iy*self.nx+ix
                self.pntx_list[ibin] = pntx
                self.pnty_list[ibin] = pnty

                
        ### initialize detector accordinate vector
        self.vdet_list = telcoord.detmm2vdet(self.pntx_list, self.pnty_list, teldef)

        
        print("initializing XRBKG")

        ### initizalize XRBKG object
        for iy in range(self.ny) :
            for ix in range(self.nx) :
                print("iy/ny = %d/%d, iy/ny = %d/%d"%(iy, self.ny, ix, self.nx)) 
                targetName = "xrbkg_x%dy%d"%(ix, iy)
                self.xrbkspec_list.append(XraySpectrum(targetName, self.caldb.rmf, self.caldb.arf))

        print("initialize XRBKG done")

        return 1



    def initXRBKSpec(self, hpspecfname) :
        if len(self.xrbkspec_list)==0 :
            if self.init_xrbkspec_list()==0 :
                ### Error case
                return 0
        
        #attrot = R.from_quat(qparam)
        vsky_list = telcoord.vdet2vsky(self.vdet_list, self.attrot, self.caldb.teldef.alignrot)
        
        ### open healpix map
        hdul = pyfits.open(hpspecfname)
        hdu = hdul[1]
        nside = hdu.header["NSIDE"]
        #hmap = hdu.data.field(icol)
        
        print("Reading spectral data from %s"%(hpspecfname))
        #nside = 64
        ipix_list = healpy.pixelfunc.vec2pix(nside, vsky_list[0], vsky_list[1], vsky_list[2])
        for iy in range(self.ny) :
            for ix in range(self.nx) :
                ibin = iy*self.nx + ix
                ipix = ipix_list[ibin]
                vSourceSpec = hdu.data[ipix][0].astype("d")
                #self.imgdat[iy][ix] = hmap[ipix]
                self.xrbkspec_list[ibin].setSourceSpec(vSourceSpec)
        hdul.close()
        return 1



    
    #def genEvents(self, tstart, tstop) :
    def simEventGen_XRBK(self, tstart, tstop) :
        if len(self.xrbkspec_list)==0 :
            print("Error in HZXSimXRBK::simEventGen_XRBK. xrbkspec_list is not initialized.")
            return 0
        
        exposure = tstop - tstart

        ### generate Events for each segments
        for iy in range(self.ny) :
            for ix in range(self.nx) :
                
                ibin = iy*self.nx + ix
                for iband in range(self.hzxpsf.nbands) :
                    elo = self.hzxpsf.elo_list[iband]
                    ehi = self.hzxpsf.ehi_list[iband]
                    
                    ncnts_expect = self.xrbkspec_list[ibin].calcObsRate(elo, ehi) * self.solidangle * exposure
                    ncnts = np.random.poisson(ncnts_expect)
                    print("ix=%d, iy=%d, expected cnts = %lf, generated cnts= %d"%(ix, iy, ncnts_expect, ncnts))
                    if ncnts<= 0 :
                        continue

                    
                    temp_vene  = self.xrbkspec_list[ibin].genObsEvents(ncnts)
                    temp_vtime = tstart + np.random.random(ncnts)*exposure

                    temp_voptx = self.pntx_list[ibin] + (np.random.random(ncnts)-0.5)*self.deltx ### shift random number range to (-0.5,0.5)
                    temp_vopty = self.pnty_list[ibin] + (np.random.random(ncnts)-0.5)*self.delty 

                    vidx =  (self.hzxpsf.nentries_list[iband] * np.random.random(ncnts)).astype('i')

                    temp_vposx = self.hzxpsf.vposx_list[iband][vidx]
                    temp_vposy = self.hzxpsf.vposy_list[iband][vidx]
                    temp_vrefx = self.hzxpsf.vrefx_list[iband][vidx]
                    temp_vrefy = self.hzxpsf.vrefy_list[iband][vidx]

                    temp_vdetx = temp_voptx + temp_vposx 
                    temp_vdety = temp_vopty + temp_vposy 
                    
                    ### Select MPO area
                    #vcut_xrtreg = np.logical_and(np.fabs(vrefx+2*optxmm)<self.xrthxy, np.fabs(vrefy+2*optymm)<self.xrthxy)
                    vcut_xrtreg = np.logical_and(np.fabs(temp_vrefx+2*temp_voptx)<self.hzxpsf.xrthxy, np.fabs(temp_vrefy+2*temp_vopty)<self.hzxpsf.xrthxy)

                    ### Exclude support structure dark lane
                    #vcut_xrtex1 = np.logical_not(np.logical_and(-self.diph-2*optxmm<vrefx, vrefx<-self.dipl-2*optxmm))
                    #vcut_xrtex2 = np.logical_not(np.logical_and( self.dipl-2*optxmm<vrefx, vrefx< self.diph-2*optxmm))
                    vcut_xrtex1 = np.logical_not(np.logical_and(-self.hzxpsf.diph-2*temp_voptx<temp_vrefx, temp_vrefx<-self.hzxpsf.dipl-2*temp_voptx))
                    vcut_xrtex2 = np.logical_not(np.logical_and( self.hzxpsf.dipl-2*temp_voptx<temp_vrefx, temp_vrefx< self.hzxpsf.diph-2*temp_voptx))

                    #vcut_xrtex3 = np.logical_not(np.logical_and(-self.diph-2*optymm<vrefy, vrefy<-self.dipl-2*optymm))
                    #vcut_xrtex4 = np.logical_not(np.logical_and( self.dipl-2*optymm<vrefy, vrefy< self.diph-2*optymm))
                    vcut_xrtex3 = np.logical_not(np.logical_and(-self.hzxpsf.diph-2*temp_vopty<temp_vrefy, temp_vrefy<-self.hzxpsf.dipl-2*temp_vopty))
                    vcut_xrtex4 = np.logical_not(np.logical_and( self.hzxpsf.dipl-2*temp_vopty<temp_vrefy, temp_vrefy< self.hzxpsf.diph-2*temp_vopty))
                    

                    #vcut_imgreg = np.logical_and(np.fabs(vposx+optxmm)<self.imghlfx, np.fabs(vposy+optymm)<self.imghlfy)
                    vcut_imgreg = np.logical_and(np.fabs(temp_vdetx)<self.hzxpsf.imghlfx, np.fabs(temp_vdety)<self.hzxpsf.imghlfy)

                    #vcut_all = np.all([vcut_xrtreg, vcut_xrtex1, vcut_xrtex2, vcut_xrtex3, vcut_xrtex4, vcut_imgreg], axis=0)
                    vmask = np.all([vcut_xrtreg, vcut_xrtex1, vcut_xrtex2, vcut_xrtex3, vcut_xrtex4, vcut_imgreg], axis=0)


                    ncnts = len(temp_vene[vmask])
                    self.vene  = np.append(self.vene,  temp_vene[vmask])
                    self.vtime = np.append(self.vtime, temp_vtime[vmask])
                    self.vdetx = np.append(self.vdetx, temp_vdetx[vmask])
                    self.vdety = np.append(self.vdety, temp_vdety[vmask])

                    vpntx = temp_voptx[vmask]
                    vpnty = temp_vopty[vmask]
                    vpnt = telcoord.detmm2vdet(vpntx, vpnty, self.caldb.teldef)
                    vra, vdec = telcoord.vdet2radec(vpnt, self.attrot, self.caldb.teldef.alignrot)
                    
                    self.vra   = np.append(self.vra,  vra)
                    self.vdec  = np.append(self.vdec, vdec)

                    ### create empty vpha array
                    self.vpha  = np.append(self.vpha, np.zeros(ncnts, 'int16'))

                    print("ix=%d, iy=%d, eband=%.1lf-%.1lf keV,  Detected cnts = %d"%(ix, iy, elo, ehi, ncnts))

        #return evts_time, evts_ene, evts_pntx, evts_pnty, evts_detx, evts_dety
        return 1

    
        
 
if __name__=="__main__" :

    from hzxpsf import HZXPSF, get_psffilelist
    from hzxcaldb import HzxCaldb
    from evtproc import EventProc
    
    hzxsim = HZXSimXRBK() ### default

    ### Event process
    evtproc = EventProc()

    ###  Observation Time and Attitude
    #qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173] ### Crab region
    qparam = [-0.82470631,  0.00189276,  0.00301748, -0.56555001]  #### Galactic center region
    tstart  = 0.0
    tstop   = 1000.0

    ### inv=True for attiude data
    ### set inv=False for the old qparam
    #hzxsim.setAttQuat(qparam, inv=False)  
    hzxsim.setAttQuat(qparam, inv=False)  

    
    ### Observation target
    #objectName = "Crab"
    #ra  = 83.633
    #dec = 22.014

    ### HZX simulate
    objectName = "XRAYSKY"


    ### init PSF database
    #psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psffits"
    psffitsdir =  "/Users/sugizaki/work/g4work/v11.3.2/myxrtg4/xg4mposim3x3/work1/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    psf_sigma = 0.2 # mmm
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)
    

    #hpspecfname = "xrbkg_rayspec.fits"
    #hzxsimdir = "/home/sugizaki/backup_d1/ep/hzxsim"
    #hpspecfname = hzxsimdir+os.sep+"xrbcat/xrbkg_rayspec.fits"
    hpspecfname = "../simdb/xrbkg_rayspec_erobin.fits"


    detid_list = range(1,17)

    for detid in detid_list :
        print("XM%02d"%(detid))

        ### init CALDB
        teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
        rmffname   = "../refdata/hzx_erosita.rmf"
        arffname   = "../refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
        vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"
        hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
        #hzxcaldb  = HzxCaldb(detid)
        hzxsim.setCaldb(hzxcaldb)

        ### init CALDB
        #hzxcaldb = Hzxcaldb(cmosid)
        #hzxsim.setCaldb(hzxcaldb)

        #hzxsim.setXRBSpec(hpspecfname)
        hzxsim.initXRBKSpec(hpspecfname)
        
        #tempfname = "temp_xrbkg_xm%02d.fits"%(detid)
        outfname = "temp_xrbkg_xm%02d.evt"%(detid)
        if os.path.isfile(outfname) :
            os.remove(outfname)
        
        hzxsim.initEventData()
        hzxsim.simEventGen_XRBK(tstart, tstop) 
        hzxsim.simDetectorResp()
        ### process
        evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, inv=False)
        evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)



