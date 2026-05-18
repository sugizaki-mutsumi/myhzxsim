#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
import matplotlib.colors as colors

### for detecor coordinate 
from teldef import Teldef
import telcoord


fig, axs = plt.subplots(3,1, figsize=(8,6), sharex=True)
plt.subplots_adjust(top=0.97, bottom=0.08, hspace=0.05)

#for ax in axs :
axs[1].set_ylabel("Counts s$^{-1}$", labelpad=20)
axs[2].set_xlabel("Time (s)")


hzxsimdir = os.getenv("HZXSIMDIR")


### set up event data
obsid = "100000"  ### Crab
#obsid = "200000"  ### test for GRBs only

teldefverid = 2 ### teldef version
dirname = obsid 
#dirname = obsid+".def"

grb_radec_detid_list = [
    [(83., 13.), 12],
    [(73., 13.), 10],
    [(83.,  3.), 11],
]


ngrbs = len(grb_radec_detid_list)
for igrb in range(ngrbs) :
    ax1 = axs[igrb]
    lcfname = "grb%d_lc.fits"%(igrb+1)

    radec = grb_radec_detid_list[igrb][0]
    detid = grb_radec_detid_list[igrb][1]
    ra_obj  = radec[0]
    dec_obj = radec[1]

    evtfname = dirname+os.sep+"hzx%sxm%02d.evt"%(obsid, detid)
    hdul = pyfits.open(evtfname)
    hdu = hdul[1]
    q1 = hdu.header["QPARAM1"]
    q2 = hdu.header["QPARAM2"]
    q3 = hdu.header["QPARAM3"]
    q4 = hdu.header["QPARAM4"]
    attrot = telcoord.attquat2rot((q1,q2,q3,q4))

    expos = hdu.header["EXPOSURE"]

    vtime = hdu.data.field("TIME") 
    vx = hdu.data.field("X") 
    vy = hdu.data.field("Y") 
    vdetx = hdu.data.field("DETX") 
    vdety = hdu.data.field("DETY") 
    vpi   = hdu.data.field("PI") 
    vra_obj  = hdu.data.field("RA_OBJ") 
    vdec_obj = hdu.data.field("DEC_OBJ") 
    vphene = hdu.data.field("PHOTON_E")

    ### cut for non x-ray bkg, soft x-ray bkg
    vcut_nxb  = vphene<0
    vcut_obj  = np.logical_and(np.fabs(vra_obj-ra_obj)<0.0001, np.fabs(vra_obj-ra_obj)<0.0001)
    vcut_sxrb = np.logical_and(np.logical_not(vcut_nxb), np.logical_not(vcut_obj))  

    ## detector region cut
    teldeffname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v%03d.fits"%(detid, teldefverid)
    this_teldef = Teldef(teldeffname) 
    detx, dety = telcoord.radec2detpix(ra_obj, dec_obj, attrot, this_teldef)
    print("Target position (DETX, DETY) = (%lf, %lf)"%(detx, dety))

    #regrad = 10.0 ### pixel (1 pixel = 0.107 mm = 1,22 arcmin 
    regrad = 20.0 ### pixel
    #regrad = 30.0 ### pixel 
    #regrad = 500.0 ### pixel, no cut
    vcut_detreg  = (vdetx-detx)**2 +(vdety-dety)**2 < regrad**2  

    vcut_nxb  = np.logical_and(vcut_nxb,  vcut_detreg)
    vcut_sxrb = np.logical_and(vcut_sxrb, vcut_detreg)

    ### light curve
    #tbinsize = 10 ## seconds
    tbinsize = 2 ## seconds
    #tbinsize = 100 ## seconds
    tstart = hdu.header["TSTART"]
    tstop  = hdu.header["TSTOP"]
    nbins = np.ceil((tstop-tstart)/tbinsize).astype(int)
    tmin = tstart
    tmax = tstart + tbinsize*nbins


    ### target object light curve
    lchist, xedges = np.histogram(vtime[vcut_detreg], nbins, (tmin, tmax))
    vx  = np.arange(nbins)*tbinsize + tbinsize/2.
    vy  = lchist/tbinsize
    vye = np.sqrt(lchist)/tbinsize
    ax1.errorbar(vx, vy, xerr=tbinsize/2., yerr=vye, fmt='.', ms=3, lw=1, label='Source events including BKG') 

    ### soft X-ray background
    lchist, xedges = np.histogram(vtime[vcut_sxrb], nbins, (tmin, tmax))
    vx  = np.arange(nbins)*tbinsize + tbinsize/2.
    vy  = lchist/tbinsize
    vye = np.sqrt(lchist)/tbinsize
    ax1.errorbar(vx, vy, xerr=tbinsize/2., yerr=vye, fmt='.', ms=3, lw=1, label='SXRB') 

    ### NXB
    lchist, xedges = np.histogram(vtime[vcut_nxb], nbins, (tmin, tmax))
    vx  = np.arange(nbins)*tbinsize + tbinsize/2.
    vy  = lchist/tbinsize
    vye = np.sqrt(lchist)/tbinsize
    ax1.errorbar(vx, vy, xerr=tbinsize/2., yerr=vye, fmt='.', ms=3, lw=1, label='NXB') 

    hdul.close()

    
    ### light curve model
    hdul = pyfits.open(lcfname)
    hdu = hdul[1]
    vlctime = hdu.data.field("TIME")
    vlcrate = hdu.data.field("RATE")
    areafact = 0.3
    ax1.plot(vlctime, vlcrate*areafact, lw=1, color='r', label='Lightcurve model') 
    hdul.close()

    ax1.text(0.02, 0.88, "GRB%d"%(igrb+1), transform=ax1.transAxes)
    if igrb==0 :
        ax1.legend()

    
axs[0].set_xlim(tmin, tmax)
    
plt.ion()
plt.show()

