#!/usr/bin/env python

import os, sys
import numpy as np
#
import astropy.io.fits as pyfits
import xspec

### CALDB
### import from $HZXSIMDIR/pylib
from respmatrix import AncillaryResponse, ResponseMatrix
from xraysource import XraySource, XraySpectrum

import matplotlib.pyplot as plt

if __name__=="__main__" :


    #fig1, ax1 = plt.subplots(figsize=(7.2, 4.8))
    fig, axs = plt.subplots(2,1, figsize=(6,6), sharex=True)
    plt.subplots_adjust(top=0.95, bottom=0.08, hspace=0.02)
    fig.align_labels()
    
    for ax in axs :
        ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlim(0.2, 10.)
        ax.tick_params(axis='both', which='both', direction='in', top=True, bottom=True, left=True, right=True)
        
    axs[0].set_ylabel(r"Photons cm$^{-1}$ s$^{-1}$ keV$^{-1}$", labelpad=12)
    axs[1].set_ylabel(r"Counts s$^{-1}$ keV$^{-1}$", labelpad=12)
    axs[1].set_xlabel("Energy (keV)")

    axs[0].set_ylim(0., 7)
    axs[1].set_ylim(0., 37)
    
    ax0 = axs[0]
    ax1 = axs[1]


    ### response files from HZXSIM
    hzxsimdir = os.getenv("HZXSIMDIR")
    RMFfname = hzxsimdir+os.sep+"refdata/hzx_erosita.rmf"

    ### Xspec Crab spectral model
    xsmodstr1 = "tbabs*power"
    xsparams1 = (0.3, 2.1, 9.7)
    objectName1 = "Object1"

    #rn = 0.0 ### surface roughness
    rn = 2.0 ### surface roughness

    region_type_list = ["imgall", "cross", "focus"]
    eband_list = [
        [0.1, 10.], [0.4, 4.0],
    ]

    for region_type in region_type_list :
    
        #ARFfname = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
        #ARFfname = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn0.0nm_cross.arf"
        #ARFfname = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn0.0nm_focus.arf"
        #ARFfname = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn2.0nm_focus.arf"
        ARFfname = hzxsimdir+os.sep+"refdata/hzxarf_withframe_rn%.1lfnm_%s.arf"%(rn, region_type)


        rmf = ResponseMatrix(RMFfname)
        arf = AncillaryResponse(ARFfname)
        srcspec  = XraySpectrum(objectName1, rmf, arf, xsmodel=xsmodstr1, xsparam=xsparams1, rmfnorm=0)

        vbinsize = srcspec.velh[1:]-srcspec.velh[0:-1]
    
        print("###  Region type", region_type)
        for eband in eband_list :
            elo = eband[0]
            ehi = eband[1]
            pflux = srcspec.calcPhotonFlux(elo, ehi)
            crate = srcspec.calcObsRate(elo, ehi)
            #print("###  Energy band %.1lf-%.1lf"%(elo, ehi))
            print("Eband %.1lf-%.1lf: Photon flux = %.4lf,  Count rate = %.4lf"%(elo, ehi, pflux, crate))
            
        #ax1.stairs(srcspec.vSourceSpec/vbinsize, srcspec.velh)
        ax1.stairs(srcspec.vFoldedSpec/vbinsize, srcspec.velh, label='%s'%(region_type))

    
    ax0.stairs(srcspec.vSourceSpec/vbinsize, srcspec.velh, color='k', label='Source spectrum')
    ax0.legend()

    ax1.legend(title='Observed spectra')
    ax1.text(0.02, 0.92, "Roughness = %.1lf nm"%(rn), transform=ax1.transAxes)

    plt.ion()
    plt.show()

    
