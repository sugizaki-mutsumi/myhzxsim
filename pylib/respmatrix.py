#!/usr/bin/env python

#import array
import os
import numpy as np
import astropy.io.fits as pyfits
#import xspec

#import datetime
#from astropy.time import Time
#import eptime


class AncillaryResponse :
    def __init__(self, arffname) :
        self.arffname = arffname
        hdul = pyfits.open(arffname)
        hdu  = hdul["SPECRESP"]
        #hdu  = hdul[1]
        self.ne   = hdu.header["NAXIS2"]
        self.vel  = hdu.data.field("ENERG_LO").astype("d")
        self.veh  = hdu.data.field("ENERG_HI").astype("d")
        self.vrsp = hdu.data.field("SPECRESP").astype("d")
        hdul.close()

        return
        

class ResponseMatrix :
    def __init__(self, rmffname) :

        self.rmffname = rmffname
        hdul = pyfits.open(rmffname)
        hdu = hdul["MATRIX"]
        detnam = hdu.header["DETNAM"]
        self.ne   = hdu.header["NAXIS2"]
        self.nchans = hdu.header["DETCHANS"]
        self.vel  = hdu.data.field("ENERG_LO").astype("d")
        self.veh  = hdu.data.field("ENERG_HI").astype("d")
        
        ### Source energy spectrum histogram
        #ne   = 1980
        self.elo = self.vel[0]
        self.ehi = self.veh[self.ne-1]

        ### Response data
        #self.matrix = hdu.data.field("MATRIX").astype('d')
        self.matrix = np.zeros((self.ne, self.nchans))

        vmat   = hdu.data.field('MATRIX')
        vfchan = hdu.data.field('F_CHAN')
        vnchan = hdu.data.field('N_CHAN')
        vngrp  = hdu.data.field('N_GRP')
        
        for irow in range(self.ne) :
            ngrp = vngrp[irow]
            for igrp in range(ngrp) :
                fchan = vfchan[irow][igrp]
                nchan = vnchan[irow][igrp]
                self.matrix[irow][fchan-1:fchan+nchan-1] = vmat[irow]

        
        ### create cdf
        self.respcdf = np.zeros((self.ne, self.nchans))
        for irow in range(self.ne) :
            self.respcdf[irow] = np.cumsum(self.matrix[irow])
            self.respcdf[irow] /= self.respcdf[irow][-1]
            

        ### Detector channel
        hdu = hdul["EBOUNDS"]
        self.vchan = hdu.data.field("CHANNEL").astype('i')
        self.vchan_emin = hdu.data.field("E_MIN").astype('d')
        self.vchan_emax = hdu.data.field("E_MAX").astype('d')
        self.chan_emin = self.vchan_emin[0] 
        self.chan_emax = self.vchan_emax[-1] 
        
        hdul.close()
        return


    def clear(self) :
        return


### move to event process code
#def procEvents(respmat, evts_ene) :
#    nevts = len(evts_ene)
#    evts_pi = np.zeros(nevts, 'int16')
#    vrand =  np.random.rand(nevts)
#    viene = ((evts_ene - respmat.elo)*respmat.ne / (respmat.ehi - respmat.elo)).astype('i')
#    np.clip(viene, 0, respmat.ne-1, None) ### cut lower/upper over limit
#    for idx in range(nevts) :
#        iene = viene[idx]
#        ichan = np.searchsorted(respmat.respcdf[iene], vrand[idx])
#        evts_pi[idx] = ichan
#    return evts_pi


if __name__=="__main__" :



    ### default RMF
    rmffname   = "../refdata/hzx_erosita.rmf"
    rmf = ResponseMatrix(rmffname)
    print("Reading response matrix done.")

    ### default ARF file
    #roughness = 0.0
    roughness = 2.0
    arffname   = "../refdata/hzxarf_withframe_rn%.1lfnm_imgall.arf"%(roughness)
    arf = AncillaryResponse(arffname)
    arffname   = "../refdata/hzxarf_withframe_rn%.1lfnm_cross.arf"%(roughness)
    arf2 = AncillaryResponse(arffname)
    arffname   = "../refdata/hzxarf_withframe_rn%.1lfnm_focus.arf"%(roughness)
    arf3 = AncillaryResponse(arffname)
    

    
    ### plot RMF data ###
    import matplotlib.pyplot as plt
    fig1, ax1 = plt.subplots(figsize=(7.2, 4.8))

    vedges = np.zeros(rmf.nchans+1)
    vedges[:-1] = rmf.vchan_emin
    vedges[-1]  = rmf.vchan_emax[-1]
    vspec = np.zeros(rmf.nchans)

    ekev_list = [0.5, 1.0, 2.0, 4.0, 5.89, 6.49, 10.0]
    for ekev in ekev_list :
        irow = np.argmin(np.fabs(rmf.vel-ekev))
        print("Photon E=%.1lf, irow=%d"%(ekev, irow))
        ax1.stairs(rmf.matrix[irow], vedges, label="%.2lf keV"%(ekev))

    ax1.set_yscale('log')
    ax1.set_ylim(1e-6, 0.2)
    ax1.set_xlim(0, 12)
    ax1.legend(fontsize=9, bbox_to_anchor=(0.9, 1), loc='upper left')
    ax1.set_xlabel("Energy (keV)")


    ### plot ARF data ###
    fig2, ax2 = plt.subplots(figsize=(7.2, 4.8))
    vedges = np.zeros(arf.ne+1)
    vedges[:-1] = arf.vel
    vedges[-1]  = arf.veh[-1]
    ax2.stairs(arf.vrsp, vedges, label="all")
    ax2.stairs(arf2.vrsp, vedges,label="corss")
    ax2.stairs(arf3.vrsp, vedges, label="focus")

    #ax1.set_yscale('log')
    #ax2.set_ylim(0, 10)
    yma = max(arf.vrsp)
    ax2.set_ylim(0.0, yma*1.1)
    ax2.set_xlim(0, 12)
    ax2.legend(fontsize=9)
    ax2.set_xlabel("Energy (keV)")
    ax2.set_ylabel(r"Effective area (cm$^{2}$)")


    plt.ion()
    plt.show()

    
