#!/usr/bin/env python

import os
#import array
import numpy as np
import astropy.io.fits as pyfits
import xspec

#import datetime
from astropy.time import Time, TimeDelta

from hzxcaldb import HzxCaldb


#define mission time default 
#from mstime import MissionTime
#MSTime = MissionTime()
#mstime_start = MSTime.asptime2mstime(Time('2025-01-01T00:00:00', scale='utc')) 
#mstime_stop  = MSTime.asptime2mstime(Time('2045-01-01T00:00:00', scale='utc')) 
#mstime_min = mstime_start
#mstime_max = mstime_stop
### mission time
from mstime import MSTime, mstime_start, mstime_stop


class XraySpectrum :
    def __init__(self, objectName, rmf, arf=None, xsmodel=None, xsparam=None, vSourceSpec=None, rmfnorm=1) :

        self.ne  = rmf.ne
        self.vel = rmf.vel
        self.veh = rmf.veh
        self.vewid = (self.veh - self.vel)
        self.vene  = (self.vel + self.veh)/2.
        if rmfnorm==1 :
            self.vrsp = np.ones(self.ne, "d")
        else :
            self.vrsp = rmf.matrix.sum(axis=1)
        rmffname = rmf.rmffname
            
        if arf!=None :
            arffname = arf.arffname
            if self.ne == arf.ne :
                self.vrsp *= arf.vrsp
            else :
                print("Warining: Energy bins of rmf and arf files do not agree.")
        else :
            arffname = ""

            
        #### set Xspec Model energy bins according to RMF
        #### This program uses only xspec.Moldel in pyXspec
        #### Responses are calcualated independently
        fs1 = xspec.FakeitSettings(rmffname, arffname, exposure=1000.)
        fs1.fileName = "temp.fak" ### file is not created. 
        if os.path.isfile(fs1.fileName) : os.remove(fs1.fileName)
        m1 = xspec.Model("power")  ### test model
        xspec.AllData.fakeit(1, fs1, noWrite=True)

        self.velh = np.array(m1.energies(1), 'd')

        self.vSourceSpec = np.zeros(self.ne)  ### create empty array
        self.vFoldedSpec = np.zeros(self.ne)  ### create empty array

        self.elo = self.vel[0]
        self.ehi = self.veh[self.ne-1]


        ### setXSpecModel
        if xsmodel!=None and xsparam!=None :
            self.setXspecModel(xsmodel, xsparam)


        elif isinstance(vSourceSpec, np.ndarray) :
            self.setSourceSpec(vSourceSpec)
        else :
            #print("Warning in XraySpectrum class. Spectrum data is not initialized.")
            pass
        
        return


    def clear(self) :
        return
    
    
    #def resetResponse(self,  rmf, arf=None, rmfnorm=1) :
    def resetSpecresp(self, rmf, arf=None, rmfnorm=1) :

        ### Reset self.vrep values
        #self.ne  = rmf.ne
        #self.vel = rmf.vel
        #self.veh = rmf.veh
        #self.vewid = (self.veh - self.vel)
        #self.vene  = (self.vel + self.veh)/2.
        if rmfnorm==1 :
            self.vrsp = np.ones(self.ne, "d")
        else :
            self.vrsp = rmf.matrix.sum(axis=1)
        
        if arf!=None :
            if self.ne == arf.ne :
                self.vrsp *= arf.vrsp
            else :
                print("Warining: Energy bins of rmf and arf files do not agree.")
        return 

        
    def setXspecModel(self, modelstr="phabs*power", params=(0.3, 2.1, 1.0) ) :

        #### Set bins according to RMF
        ### Model spectrum
        m1 = xspec.Model(modelstr)
        m1.setPars(params)
        m1.show()
        
        self.vSourceSpec = np.array(m1.values(1), 'd')
        self.vFoldedSpec = self.vSourceSpec * self.vrsp

        
    def setSourceSpec(self, vSourceSpec) :
        self.vSourceSpec = vSourceSpec
        self.vFoldedSpec = self.vSourceSpec * self.vrsp
        return


    def getSourceSpec(self) :
        return self.vSourceSpec

        
    def setPhabsPowerlaw(self, nH, phoIndex, norm=1.0) :
        xsmodel = "phabs*power"
        self.setXspecModel(xsmodel, (nH, phoIndex, norm) )
        return
    

    def genSourceEvents(self, num, elo=0.1, ehi=10.0) :

        ### cut lower/higher energy band
        this_spec = self.vSourceSpec * ((elo<=self.vel)*(self.veh<=veh))
        
        cdf = np.cumsum(this_spec)    # 
        cdf = cdf / cdf[-1]       # cdf0 value from h[0] to 1.0
        cdf0 = np.append(0.0, cdf[:-1])  ### index shift -1.  cdf0 value from 0.0 to 1.0 - h[-1]
        cdfdel = np.diff(cdf, prepend=0.0) ### nelem = len(xbins)-1

        vrand = np.random.rand(num)
        ibins = np.searchsorted(cdf, vrand, 'right') 
        vrandel = (vrand-cdf0[ibins])/(cdfdel[ibins])
        values = self.vel[ibins] + vrandel*self.vewid[ibins]

        return values
        
        

    def genObsEvents(self, num, elo=0.1, ehi=10.0, vfact=1.0) :

        if isinstance(vfact, np.ndarray) :
            is_vfact = True
        else :
            is_vfact = False
            
        if is_vfact and len(vfact)!=self.ne :
            print("Error in XraySpectrum::genObsEvents. len(vfact)!= Num of energy bins.")
            return

        ### cut lower/higher energy band
        #this_spec = self.vFoldedSpec * vfact * ((elo<=self.vel)*(self.veh<=ehi))
        ### need check later
        this_spec = self.vFoldedSpec * vfact * ((elo<=self.vel)*(self.vel<ehi))
        
        cdf = np.cumsum(this_spec)    # 
        cdf = cdf / cdf[-1]       # cdf0 value from h[0] to 1.0
        cdf0 = np.append(0.0, cdf[:-1])  ### index shift -1.  cdf0 value from 0.0 to 1.0 - h[-1]
        cdfdel = np.diff(cdf, prepend=0.0) ### nelem = len(xbins)-1

        vrand = np.random.rand(num)
        ibins = np.searchsorted(cdf, vrand, 'right') 
        vrandel = (vrand-cdf0[ibins])/(cdfdel[ibins])

        values = self.vel[ibins] + vrandel*self.vewid[ibins]

        return values


    def calcPhotonFlux(self, elo=0.1, ehi=10.0) :
        pflux = (self.vSourceSpec * ((elo<=self.vel)*(self.veh<=ehi))).sum()
        return pflux

        
    def calcObsRate(self, elo=0.1, ehi=10.0, vfact=1.0) :
        crate = (self.vFoldedSpec * vfact * ((elo<=self.vel)*(self.veh<=ehi))).sum()
        return crate
        

class LightCurve :
    def __init__(self, objectName, obsRate=0.0, tstart=mstime_start, tstop=mstime_stop) :

        self.objectName = objectName
        self.tstart = tstart
        self.tstop  = tstop

        ### case 1: constant flux 
        self.constObsRate = obsRate

        ### case 2: data with 1-D array
        self.vtime = None
        self.vrate = None
        
        ### case 3: model function 
        #self.lcfunc = None

        return


    def clear(self) :
        pass
        
    def setConstObsRate(self, obsRate) :
        self.constObsRate = obsRate

        self.vtime = None
        self.vrate = None

    def setLcdata(self, vtime, vrate) :
        ### vrate means count rate in 0.1-10 keV band

        self.tstart = vtime[0]
        self.tstop  = vtime[-1]
        self.vtime = vtime
        self.vrate = vrate

        self.constObsRate = None
        return

    
    #def readLcfile(self, lcdatfname, hdunum=1, timecol="MJD", ratecol="RATE", convfact=1.0) :
    def readLcfile(self, lcdatfname, hdunum, timecol, ratecol, convfact) :
        ### read LC fits file
        ### format talbe
        hdul = pyfits.open(lcdatfname)
        hdu = hdul[hdunum]
        if not (timecol in hdu.columns.names and ratecol in hdu.columns.names ) :
            print("lcfname=%s timecol=%s ratecol=%s do not exist"%(lcdatfname, timecol, ratecol))
            return
        
        if timecol=="MJD" :          
            self.vtime = MSTime.mjdtt2mstime(hdu.data.field(timecol).astype('d'))
        else :
            self.vtime = hdu.data.field(timecol).astype('d')

        self.vrate = hdu.data.field(ratecol).astype('d') * convfact
        hdul.close()

        self.tstart = self.vtime[0]
        self.tstop  = self.vtime[-1]


        ###  Clear other LC parameters
        self.constObsRate = None

        
    """
    def setLcfunc(self, lcfunc, params, tstart, tstop) :
        self.tstart = tstart
        self.tstop  = tstop
        self.lcfunc = ROOT.TF1("lcfunc_%s"%(self.objectName), lcfunc, tstart, tstop)
        nparams = len(params)
        for idx in range(nparams) :
            self.lcfunc.SetParameter(idx, params[idx])
    
        self.constObsRate = None
        self.vtime   = None
        self.vrate   = None
    """

    
    def getObsRate(self, evtime) :
        #print("evtime = %lf, tstart = %lf, tstop= %lf"%(evtime, self.tstart, self.tstop))
        obsRate = 0.0
        if self.tstart<=evtime and evtime<=self.tstop :
            if self.constObsRate != None :
                obsRate = self.constObsRate
            #elif self.lcfunc != None :
            #    obsRate = self.lcfunc.Eval(evtime)

            elif self.vtime!=None and self.vrate!=None :
                obsRate = np.interp(evtime, self.vtime, self.vrate)
                
        return obsRate
        
        
    def genObsEvents(self, tstart, tstop, fact=1.0) :

        #print(tstart, tstop, self.tstart, self.tstop)
        if self.tstop<=tstart or tstop<=self.tstart :
            return np.zeros(0)
            
        obs_tstart = max(self.tstart, tstart)
        obs_tstop  = min(self.tstop, tstop)
        exposure = obs_tstop - obs_tstart
        if exposure <=0.0 :
            return np.zeros(0)

        #print(obs_tstart, obs_tstop, exposure)

        if self.constObsRate != None :
            ncnts_expect = self.constObsRate*exposure
            ncnts = np.random.poisson(ncnts_expect * fact)
            evts_time = np.random.random(ncnts)*exposure + obs_tstart

        
        #elif self.vtime!=None and self.vrate!=None :
        elif isinstance(self.vtime, np.ndarray) and isinstance(self.vrate, np.ndarray) :

            linear_approximation = True
            #linear_approximation = False
            
            if not linear_approximation :
                ### Step-bin approximation
                nbin = int(exposure) ## binsize ~ 1 sec
                tedges = np.linspace(obs_tstart, obs_tstop, nbin+1)
                vtdiff = np.diff(tedges)
                vtime  = tedges[:-1] + vtdiff/2.0
                vrate  = np.interp(vtime, self.vtime, self.vrate)  
                ncnts_expect  = (vrate*vtdiff).sum()
                if ncnts_expect<=0.0 :
                    return np.zeros(0)
            
                ncnts = np.random.poisson(ncnts_expect * fact)
                if ncnts==0 :
                    return np.zeros(0)
            
                evts_time = np.zeros(ncnts)
                ### filter negative value in vrate
                vrate *= (vrate>0.0)
            
                cdf = np.cumsum(vrate)
                cdf = cdf / cdf[-1]
                cdf0 = np.append(0.0, cdf[:-1])
                cdfdel = np.diff(cdf, prepend=0.0)
                
                vrand = np.random.rand(ncnts)
                ibins = np.searchsorted(cdf, vrand, 'right')
                vrandel = (vrand-cdf0[ibins])/(cdfdel[ibins])
                
                evts_time = tedges[ibins] + vrandel*vtdiff[ibins]
            
            else :
                ### Linear-function approximation
                vidx_gti = np.logical_and(obs_tstart<=self.vtime, self.vtime<=obs_tstop)
                if not np.any(vidx_gti) :
                    rate_tstart = np.interp(obs_tstart, self.vtime, self.vrate)
                    rate_tstop  = np.interp(obs_tstop, self.vtime, self.vrate)
                    vtime = np.array((obs_tstart,obs_tstop), 'd')
                    vrate = np.array((rate_tstart,rate_tstop), 'd')
                    
                else :
                    vtime = self.vtime[vidx_gti]
                    vrate = self.vrate[vidx_gti]
                    if obs_tstart<vtime[0] :
                        rate_tstart = np.interp(obs_tstart, self.vtime, self.vrate)
                        vtime = np.append(obs_tstart, vtime)
                        vrate = np.append(rate_tstart, vrate)
                    if vtime[-1]< obs_tstop :
                        rate_tstop = np.interp(obs_tstop, self.vtime, self.vrate)
                        vtime = np.append(vtime, obs_tstop)
                        vrate = np.append(vrate, rate_tstop)

                vx = vtime
                vy = vrate
                vdx = np.diff(vx) ### ndim = nv-1
                vdy = np.diff(vy)

                ### sum of (event rate * dt)
                ncnts_expect  = ((vy[1:]+vy[:-1])*vdx).sum()/2.0
                if ncnts_expect<=0.0 :
                    return np.zeros(0)
            
                ncnts = np.random.poisson(ncnts_expect * fact)
                if ncnts==0 :
                    return np.zeros(0)

                ### filter negative value in vrate
                vy *= (vy>0.0)

                va = vdy/vdx
                vb = np.copy(vy[:-1])  ### ndim = nv-1

                cdf = np.cumsum(0.5*va*vdx*vdx + vb*vdx)
                norm = cdf[-1]
                
                #print("normalization = ", norm)
                #cdf = cdf/norm                ### cdf0 value from h[0] to 1.0  ndim = nv-1
                cdf/=norm                ### cdf0 value from h[0] to 1.0  ndim = nv-1
                cdf0 = np.append(0.0, cdf[:-1])    ### index shift -1.  cdf0 value from 0.0 to 1.0 - h[-1]
                cdfdel = np.diff(cdf, prepend=0.0) ### ndim = nv-1

                va/=norm
                vb/=norm

                evtrand = np.random.rand(ncnts)
                evtvidx = np.searchsorted(cdf, evtrand, 'right') 
                evtrr = (evtrand-cdf0[evtvidx])  #/(cdfdel[ibins])
                evta = va[evtvidx] 
                evtb = vb[evtvidx] 
                evtx = vx[evtvidx]
    
                ### random number fraction
                #print(vrr)

                ### Case of a!=0.0 
                evta_zero_case = (evta==0.0)
                evta_nonzero = evta + (evta_zero_case)*1.
                xx_case1 = (-evtb + np.sqrt(evtb*evtb + 2*evta*evtrr))/evta_nonzero

                ### Case of a==0.0 
                evtb_zero_case = (evtb==0.0)
                evtb_nonzero = evtb + (evtb_zero_case)*1.
                xx_case2 = evtrr/evtb_nonzero
    
                #xx = xx_case1*(~evta_zero_case) + xx_case2*(evta_zero_case)
                values = evtx + xx_case1*(~evta_zero_case) + xx_case2*(evta_zero_case)
                evts_time = values
                #return values, norm


        evts_time.sort()
        return evts_time
            

    
class XraySource :
    def __init__(self, objectName, ra, dec, tstart, tstop, rmf, arf, xsmodel="phabs*power", xsparam=(0.3, 2.1, 1.0)) :

        self.objectName = objectName
        self.ra  = ra
        self.dec = dec
        
        ### create default light curve
        obsRate = 0.0  # default
        self.lcurve = LightCurve(objectName, obsRate, tstart, tstop) 

        ### set XraySpectrum and modify LightCurve rate
        self.spectrum = XraySpectrum(objectName, rmf, arf, xsmodel, xsparam) 
        obsRate = self.spectrum.calcObsRate()
        self.lcurve.setConstObsRate(obsRate)

    def clear(self) :
        self.spectrum.clear()
        self.lcurve.clear()
        
    def resetSpecresp(self, rmf, arf=None, rmfnorm=1) :
        self.spectrum.resetSpecresp(rmf, arf, rmfnorm)
        
    def setObjectNameRADEC(self, objectName, ra, dec) :
        self.objectName = objectName
        self.ra  = ra
        self.dec = dec
    
    def setXspecModel(self, xsmodel="phabs*power", xsparam=(0.3, 2.1, 1.0) ) :
        self.spectrum.setXspecModel(xsmodel, xsparam)
        obsRate = self.spectrum.calcObsRate()
        #print("obsRate = %lf"%(obsRate))
        self.lcurve.setConstObsRate(obsRate)

        
    #def setLcdata(self, vtime, vrate, elo, elhi) :
    def setLcdata(self, vtime, vpflux, elo, elhi) :

        ### Normalize light curve data in units of expected count rate from units of photon flux 
        pflux = self.spectrum.calcPhotonFlux(elo, ehi)
        crate = self.spectrum.calcObsRate()
        eafact = crate/pflux

        print("Photon flux -> Count rate conversion factor = %lf"%(eafact))

        self.lcurve.setLcdata(vtime, vpflux*eafact)
        return

    
    def readLcfile(self, lcdatfname, hdunum=1, timecol="TIME", ratecol="RATE", convfact=1.0, elo=0.1, ehi=10.0) :
        if convfact>0.0 :
            eafact = convfact
        else :
            ### Light curve data is in units of photon flux 
            pflux = self.spectrum.calcPhotonFlux(elo, ehi)
            crate = self.spectrum.calcObsRate()
            eafact = crate/pflux

        print("Photon flux -> Count rate conversion factor = %lf"%(eafact))

        self.lcurve.readLcfile(lcdatfname, hdunum, timecol, ratecol, eafact)
        return

    
    def getObsRate(self, evtime) :
        obsRate = self.lcurve.getObsRate(evtime)
        return obsRate
    

    #def genEvents(self, tstart, tstop, emin=0.1, emax=10.0) :
    def genEvents(self, tstart, tstop, elo=0.1, ehi=10.0, vigfact=1.0) :

        if elo<=self.spectrum.elo and self.spectrum.ehi<=ehi and not isinstance(vigfact, np.ndarray) :
            specfact = vigfact
        else :
            obsrate_total = self.spectrum.calcObsRate(self.spectrum.elo, self.spectrum.ehi, 1.0)
            obsrate_eband = self.spectrum.calcObsRate(elo, ehi, vigfact)
            specfact = obsrate_eband/obsrate_total
            print("elo=%.2lf ehi=%.2lf specfact=%lf"%(elo, ehi, specfact))
            #print("elo=%.2lf ehi=%.2lf specfact=%lf, areafact=%lf"%(elo, ehi, specfact, fact))


        
        evts_time = self.lcurve.genObsEvents(tstart, tstop, specfact)
        nevts = len(evts_time)
        print("Energy band %.2lf - %.2lf keV nevents = %d"%(elo, ehi, nevts))
        
        evts_ene = self.spectrum.genObsEvents(nevts, elo, ehi, vigfact)
        return evts_time, evts_ene
    
    

if __name__=="__main__" :

    import matplotlib.pyplot as plt

    detid = 1
    detnam = "XM%02d"%(detid)

    teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
    rmffname   = "../refdata/hzx_erosita.rmf"
    #arffname   = "../refdata/hizx_arf_all.fits"
    arffname   = "../refdata/hzxarf_withframe_rn0.0nm_focus.arf"
    vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"
    
    hzxcaldb =  HzxCaldb(1, teldefname, rmffname, arffname, vigfname)

    ### plot
    fig1, ax1 = plt.subplots(figsize=(7.2, 4.8))
    srcspec = XraySpectrum("Object1", hzxcaldb.rmf, hzxcaldb.arf, xsmodel="tbabs*power", xsparam=(0.3, 2.1, 10.0))
    #srcspec = XraySpectrum("Object1", hzxcaldb.rmf, hzxcaldb.arf, xsmodel="tbabs*power", xsparam=(0.0, 1.0, 10))

    #geoarea = 6.144*6.144 # cm^2
    vbinsize = srcspec.velh[1:]-srcspec.velh[0:-1]
    ax1.stairs(srcspec.vSourceSpec/vbinsize, srcspec.velh)
    ax1.stairs(srcspec.vFoldedSpec/vbinsize, srcspec.velh)

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.2, 12)
    #ax1.set_ylim(0.5, 50)
    ax1.set_ylim(0.02, 20)

    ax1.set_xlabel("Energy (keV)")
    ax1.set_ylabel(r"Counts s$^{-1}$ keV$^{-1}$")
    
    pflux = srcspec.calcPhotonFlux()
    crate = srcspec.calcObsRate()
    print("Photon Flux = ", pflux)
    print("Obs rate = ", crate)
    

    eband_list = [
        [0.5, 1.0],
        [1.0, 2.0],
        [2.0, 4.0],
    ]
    for eband in eband_list :
        elo = eband[0]
        ehi = eband[1]
        pflux = srcspec.calcPhotonFlux(elo, ehi)
        crate = srcspec.calcObsRate(elo, ehi)
        print("%.1lf - %.1lf keV, PhotonFlux = %lf, ObsRate = %lf"%(elo, ehi, pflux, crate))
    

    mjdtt_start = Time('2026-04-01T00:00:00', scale='tt').tt.mjd
    mjdtt_stop  = Time('2026-04-01T00:20:00', scale='tt').tt.mjd  ### 20*60 = 1200 sec
    #tstart = MSTime.mjdtt2mstime(mjdtt_start)
    #tstop  = MSTime.mjdtt2mstime(mjdtt_stop)
    tstart = 0.
    tstop  = 5000.
    expos = tstop - tstart
    print("Exposure = %lf sec.", expos)
    ra  = 83.633
    dec = 22.014
    xraysrc = XraySource("Source1", ra, dec, tstart, tstop, hzxcaldb.rmf, hzxcaldb.arf) 

    xraysrc.setXspecModel("tbabs*power", (0.3, 2.1, 10.0))
    obsRate = xraysrc.getObsRate(tstart)
    print("obsRate = %lf c/s"%(obsRate))

    evts_time, evts_ene = xraysrc.genEvents(tstart, tstop)

    #fig2, ax2 = plt.subplots(figsize=(7.2, 4.8))

    #vsimspec, xedges = np.histogram(evts_ene, srcspec.velh)
    #vx  = (xedges[:-1] + xedges[1:])/2.
    #vxe = (xedges[1:] - xedges[:-1])/2.
    #vy  = vsimspec/vbinsize/expos
    #vye = np.sqrt(vsimspec)/vbinsize/expos
    #ax2.errorbar(vx, vy, xerr=vxe, yerr=vye, fmt='.', ms=3, lw=1)

    ### rebin
    #xbins = np.logspace(np.log10(0.2),np.log10(12.), 100)
    #xbins = np.logspace(np.log10(0.2),np.log10(12.), 200)
    xbins = np.logspace(np.log10(0.2),np.log10(12.), 500)
    #xbins = np.append(xraysrc.spectrum.vel, xraysrc.spectrum.veh[-1]) 
    #xbins = xbins[::2]
    vspechist, xedges = np.histogram(evts_ene, xbins)
    vxcen = (xbins[1:] + xbins[:-1])/2
    vxwid =  xbins[1:] - xbins[:-1]
    vy  = vspechist/vxwid/expos
    vye = np.sqrt(vspechist)/vxwid/expos
    
    ax1.errorbar(vxcen, vy, xerr=vxwid/2, yerr=vye, fmt='.', ms=3, lw=1, color='r')
    
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    #ax2.set_xlim(0.2, 12)
    ##ax2.set_ylim(0.5, 50)
    #ax2.set_ylim(0.02, 20)
    
    
    plt.ion()
    plt.show()
    
