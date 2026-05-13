#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits


class HZXPSF :
    def __init__(self, eband_le=0.1, eband_he=12.0, sigma=0.2, elist=[], flist=[], psfdir="") :

        self.version_number = 1
        
        if not (len(elist)>0 and len(elist)==len(flist) ) :
            print("Error in HZXPSF. Energy band list or PSF file list are illegal!") 
            return
        
        self.elist   = np.array(elist)
        self.nbands = len(self.elist)
        self.eband_le = eband_le
        self.eband_he = eband_he
        self.elo_list = np.zeros(self.nbands)
        self.ehi_list = np.zeros(self.nbands)
        for iband in range(self.nbands) :
            if iband == 0 :
                self.elo_list[iband] = self.eband_le
            else :
                self.elo_list[iband] = (self.elist[iband-1] + self.elist[iband])/2.0

            if iband == self.nbands -1 :
                self.ehi_list[iband] = self.eband_he
            else :
                self.ehi_list[iband] = (self.elist[iband] + self.elist[iband+1])/2.0
        


        ### MPO array dimension
        ### support structure dip width, low, high
        self.dipw = 3.0 # mm
        self.dipl = 20.0  ## mm
        self.diph = self.dipl +  self.dipw
        self.xrthxy = 60.0 + self.dipw

        
        ### imager pixel dimension
        ### pixel size = 108 um
        ### 108 um * 512 = 55.296 mm ~ 55 mm
        self.pixsize = 0.108 ## mm
        self.nxpix  = 512  
        self.nypix  = 512  
        self.imghlfx = self.pixsize * self.nxpix / 2.0
        self.imghlfy = self.pixsize * self.nypix / 2.0

        ## EP WXT case
        #self.pixsize = 0.015 ## mm
        #self.imghlfx = self.pixsize*2048
        #self.imghlfy = self.pixsize*2048

        ### PSF smoothing parameter
        #sigma = 0.2  # 0.2 mm = 4 arcmin  


        #self.tree_list = []
        self.vposx_list = []
        self.vposy_list = []
        self.vrefx_list = []
        self.vrefy_list = []

        self.vpsfx_list = []
        self.vpsfy_list = []
        
        self.nentries_list = np.zeros(self.nbands, 'int32')

        #for photon_ene in self.elist :            
        for iband in range(self.nbands) :

            photon_ene = self.elist[iband]
            psffname = psfdir+os.sep+flist[iband]

            hdul = pyfits.open(psffname)
            hdu = hdul[1]
            nrows = hdu.header["NAXIS2"]

            self.nentries_list[iband] = nrows
            print("Num of total events of photon E=%lf keV = %d"%(photon_ene, self.nentries_list[iband]))
            #self.vposx_list.append(hdu.data.field("posX").astype('d'))
            #self.vposy_list.append(hdu.data.field("posY").astype('d'))
            self.vposx_list.append(hdu.data.field("posX").astype('d') + np.random.normal(0.0, sigma, nrows))
            self.vposy_list.append(hdu.data.field("posY").astype('d') + np.random.normal(0.0, sigma, nrows))
            self.vrefx_list.append(hdu.data.field("refX").astype('d'))
            self.vrefy_list.append(hdu.data.field("refY").astype('d'))
            
            hdul.close()

            self.vpsfx_list.append(None) ### create dummy content
            self.vpsfy_list.append(None)
            

        ### default parameters
        self.optxmm = -999. 
        self.optymm = -999.
        self.psfnevts_list = np.zeros(self.nbands, 'int32') ### create defualt list
        self.refnevts_list = np.zeros(self.nbands, 'int32') ### create defualt list
        self.areafact_list = np.zeros(self.nbands, 'd')     ### 64bit double
        
        self.initPSFevts(0.0, 0.0) ## initialize reference
        for iband in range(self.nbands) :
            self.refnevts_list[iband] = self.psfnevts_list[iband]


        return


    def getEneBandIdx(self, ene) :
        iband = np.argmin(np.fabs(self.elist - ene))
        return iband
        

    def initPSFevts(self, optxmm, optymm) :

        if self.optxmm == optxmm and  self.optymm == optymm :
            #return self.psfnevts
            return
        
        print("Initialize PSF for optxy = %lf, %lf"%(optxmm, optymm))


        for iband in range(self.nbands) :

            vposx = self.vposx_list[iband]
            vposy = self.vposy_list[iband]
            vrefx = self.vrefx_list[iband]
            vrefy = self.vrefy_list[iband]
            
            ### Select MPO area
            vcut_xrtreg = np.logical_and(np.fabs(vrefx+2*optxmm)<self.xrthxy, np.fabs(vrefy+2*optymm)<self.xrthxy)

            ### Exclude support structure dark lane
            vcut_xrtex1 = np.logical_not(np.logical_and(-self.diph-2*optxmm<vrefx, vrefx<-self.dipl-2*optxmm))
            vcut_xrtex2 = np.logical_not(np.logical_and( self.dipl-2*optxmm<vrefx, vrefx< self.diph-2*optxmm))

            vcut_xrtex3 = np.logical_not(np.logical_and(-self.diph-2*optymm<vrefy, vrefy<-self.dipl-2*optymm))
            vcut_xrtex4 = np.logical_not(np.logical_and( self.dipl-2*optymm<vrefy, vrefy< self.diph-2*optymm))
        
            ### select image area
            vcut_imgreg = np.logical_and(np.fabs(vposx+optxmm)<self.imghlfx, np.fabs(vposy+optymm)<self.imghlfy)
            vcut_all = np.all([vcut_xrtreg, vcut_xrtex1, vcut_xrtex2, vcut_xrtex3, vcut_xrtex4, vcut_imgreg], axis=0)
            
            self.vpsfx_list[iband] = vposx[vcut_all] + optxmm
            self.vpsfy_list[iband] = vposy[vcut_all] + optymm
            self.psfnevts_list[iband] = len(self.vpsfx_list[iband])
            if self.refnevts_list[iband]>0 :
                self.areafact_list[iband] = self.psfnevts_list[iband]/self.refnevts_list[iband]
            else :
                self.areafact_list[iband] = 1.0
                
            print("Energy band %5.2lf-%5.2lf keV nevents=%8d, areafact=%lf"%(self.elo_list[iband], self.ehi_list[iband],
                                                                             self.psfnevts_list[iband], self.areafact_list[iband]))

        self.optxmm = optxmm
        self.optymm = optymm 

        return 

    
    def countPSFevts(self, optxmm, optymm, iband=0, widx=0.6, widy=0.6, focus=True) :
        ncnts = 0
        if not (0<=iband and iband<=self.nbands) :
            return ncnts
        
        if focus :
            vcut = np.logical_and(np.fabs(self.vpsfx_list[iband]-optxmm)<widx, np.fabs(self.vpsfy_list[iband]-optymm)<widy)
        else :
            vcut = np.logical_or(np.fabs(self.vpsfx_list[iband]-optxmm)<widx, np.fabs(self.vpsfy_list[iband]-optymm)<widy)
        ncnts = sum(vcut*1)
        return ncnts
    

    def genPSFevts(self, vene) :
        num = len(vene)
        evts_vx = np.zeros(num)
        evts_vy = np.zeros(num)
        for irow in range(num) :
            ene = vene[irow] 
            iband = self.getEneBandIdx(ene)  
            iv = int(np.random.random()*self.psfnevts_list[iband])
            ###print("ene=%lf, iband=%d, iv=%d"%(ene, iband, iv))
            evts_vx[irow] = self.vpsfx_list[iband][iv]
            evts_vy[irow] = self.vpsfy_list[iband][iv]
        return evts_vx, evts_vy

    
    def genPSFevts_defframe(self, vene, detx=0.0, dety=0.0, detxmi=-30, detxma=30, detymi=-30, detyma=30) :
        self.initPSFevts(0.0, 0.0) ## initialize PSF at the default referene
        num = len(vene)
        evts_vx = np.zeros(num)
        evts_vy = np.zeros(num)
        nmax = 1000000 ### limit of maximum events generated 
        irow = 0
        #while idx < nmax :
        for idx in range(nmax) :
            ene = vene[irow] 
            iband = self.getEneBandIdx(ene)  
            iv = int(np.random.random()*self.psfnevts_list[iband])

            evtx = self.vpsfx_list[iband][iv] + detx
            evty = self.vpsfy_list[iband][iv] + dety

            ##if detxmi<=evtx and detx<=detxma and detymi<=evty and dety<=detyma :
            if detxmi<evtx and evtx<detxma and detymi<evty and evty<detyma :
                evts_vx[irow] = evtx
                evts_vy[irow] = evty
                irow +=1
            if num<=irow :
                break
            
        return evts_vx, evts_vy, irow

    
    def genPSFevts_iband(self, iband, num) :
        evts_vx = np.zeros(num)
        evts_vy = np.zeros(num)
        for irow in range(num) :
            iv = int(np.random.random()*self.psfnevts_list[iband])
            evts_vx[irow] = self.vpsfx_list[iband][iv] 
            evts_vy[irow] = self.vpsfy_list[iband][iv] 
        return evts_vx, evts_vy

    def close(self) :
        pass


    
def get_psffilelist(enelist, dirname) :
    this_elist = []
    psffilelist  = []
    nevtstr = "n1e6"
    for ene in enelist :
        fname = "gamma%.0lfev_x0y0_s60mm_fw0.0_box_rn0.0_ma0.0_%s.fits"%(ene*1e3, nevtstr)
        if os.path.isfile(dirname+os.sep+fname) :
            this_elist.append(ene)
            psffilelist.append(fname)
        
    return this_elist, psffilelist


if __name__=="__main__" :

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    #photon_ene = 0.2 ## keV
    #photon_ene = 0.5 ## keV
    photon_ene = 1.0 ## keV
    #photon_ene = 2.0 ## keV
    #photon_ene = 5.0 ## keV

    psf_sigma = 0.2 # mmm
    #psf_sigma = 0.0 # mmm
    

    #psffitsdir =  os.getenv("HZXSIMDIR")+ose.sep+"refdata/psfdb"
    psffitsdir =  "../refdata/psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )

    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)
    
    hzxpsf = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    #print(hzxpsf)
    
    optxy_list = [
        [0.0, 0.0],
        #[3.0, 0.0],
        #[5.0, 0.0],
        #[8.0, 2.0],
        [11.0, 20.0],
        [15.0, 34.0],
        #[20.0, 38.0],
        ##[0.0, 3.0],
    ]


    #nx  = 1024
    #ny  = 1024
    pixsize = hzxpsf.pixsize
    nx  = hzxpsf.nxpix  
    ny  = hzxpsf.nypix  
    xmi = -hzxpsf.imghlfx  
    xma =  hzxpsf.imghlfx  
    ymi = -hzxpsf.imghlfy  
    yma =  hzxpsf.imghlfy  

    
    #himg = ROOT.TH2F("himg", "", nx, xmi, xma, ny, ymi, yma)
    #fig,  ax  = plt.subplots(figsize=(8,4))    
    #fig2 = plt.figure(figsize=(6,6)) # , layout="constrained")
    fig, ax = plt.subplots(figsize=(5,5))
    ax.axis('equal')
    
    nphotons = 100000
    #nphotons = 1000000
    for optxy in optxy_list :
        optxmm = optxy[0]
        optymm = optxy[1]
        
        hzxpsf.initPSFevts(optxmm, optymm)
        #print("optxy=(%lf,%lf) areafact=%lf"%(optxmm, optymm, hzxpsf.areafact))  

        #for idx in range(hzxpsf.psfnevts) :
        #    himg.Fill(hzxpsf.psfv1[idx], hzxpsf.psfv2[idx]) 

        #nevts = int(nphotons * hzxpsf.areafact)
        #nevts = rand.Poisson(nphotons * hzxpsf.areafact)

        iband = hzxpsf.getEneBandIdx(photon_ene)
        areafact = hzxpsf.areafact_list[iband]
        nevts = np.random.poisson(nphotons * areafact)

        evts_optx, evts_opty = hzxpsf.genPSFevts_iband(iband, nevts)
        #for idx in range(nevts) :
        #    himg.Fill(evts_optx[idx], evts_opty[idx]) 

        hist2d, xedges, yedges = np.histogram2d(evts_optx, evts_opty, bins=(nx, ny), range=((xmi,xma), (ymi,yma)))
        hist2d = hist2d.T
        im = ax.imshow(hist2d, interpolation='nearest', origin='lower',
                       norm=colors.LogNorm(),
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])


    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    plt.ion()
    plt.show()


