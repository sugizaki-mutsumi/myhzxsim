#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
import matplotlib.colors as colors

### for detecor coordinate 
from teldef import Teldef
import telcoord


hzxsimdir = os.getenv("HZXSIMDIR")

#fig1, ax1 = plt.subplots(figsize=(6,6))
#ax2.set_yscale('log')
#ax2.set_xscale('log')

fig, axs = plt.subplots(2,2, figsize=(8,8))
#fig.subplots_adjust(left=0.05, right=0.98, bottom=0.12, top=0.92, wspace=0.25)
fig.subplots_adjust(right=0.93, hspace=0.25, wspace=0.3)



### Crab coordinates
ra_crab  = 83.63321 
dec_crab = 22.01446

### set up event data
obsid = "100000"  ### Crab
detid = 7

## teldef
teldefverid = 2 ### teldef version
teldeffname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v%03d.fits"%(detid, teldefverid)
this_teldef = Teldef(teldeffname) 



dirname_list = [
    obsid+".1st",  ### default
    obsid  
]
psfsigma_list = [
    0.2,
    0.4
]
#for dirname in dirname_list :
for idx in range(2) :
    dirname  = dirname_list[idx] 
    psfsigma = psfsigma_list[idx] 

    evtfname = dirname+os.sep+"hzx%sxm%02d.evt"%(obsid, detid)
    hdul = pyfits.open(evtfname)
    hdu = hdul[1]
    #expos = 1000.
    q1 = hdu.header["QPARAM1"]
    q2 = hdu.header["QPARAM2"]
    q3 = hdu.header["QPARAM3"]
    q4 = hdu.header["QPARAM4"]
    attrot = telcoord.attquat2rot((q1,q2,q3,q4))
    
    expos = hdu.header["EXPOSURE"]
    vx = hdu.data.field("X") 
    vy = hdu.data.field("Y") 
    vdetx = hdu.data.field("DETX") 
    vdety = hdu.data.field("DETY") 
    vpi   = hdu.data.field("PI") 
    vra_obj  = hdu.data.field("RA_OBJ") 
    vdec_obj = hdu.data.field("DEC_OBJ") 
    vphene = hdu.data.field("PHOTON_E")
    hdul.close()
    

    ### get DETXY position
    detx_crab, dety_crab = telcoord.radec2detpix(ra_crab, dec_crab, attrot, this_teldef)
    print("Crab position (DETX, DETY) = (%lf, %lf)"%(detx_crab, dety_crab))
    
    
    #### Crab cut
    vcut_crab = np.logical_and(np.fabs(vra_obj-ra_crab)<0.0001, np.fabs(vra_obj-ra_crab)<0.0001)
    
    ### select for nxb, x-ray backround
    vcut_nxb  = vphene<0
    vcut_sxrb = np.logical_and(np.logical_not(vcut_nxb), np.logical_not(vcut_crab))  
    
    
    #nx = this_teldef.det_xsiz
    detxmi = this_teldef.detxpix_min
    detxma = this_teldef.detxpix_max
    detymi = this_teldef.detypix_min
    detyma = this_teldef.detypix_max
    #detxhist, xedges = np.histogram(vdetx, nx, (xmi, xma))


    vdelx = vdetx-detx_crab
    vdely = vdety-dety_crab
    #delxysiz = 101
    delxysiz = 51
    #delxymi = -delxysiz/2.0 
    #delxyma = +delxysiz/2.0 
    #delxybins = np.linspace(delxymi, delxyma, delxysiz+1)
    vcut = np.all([np.fabs(vdelx)<delxysiz/2.0, np.fabs(vdely)<delxysiz/2.0], axis=0)
    

    axs[0][idx].plot(vdetx, vdety, '.', ms=0.1, color='k')
    axs[0][idx].set_xlim(detxmi, detxma)
    axs[0][idx].set_ylim(detymi, detyma)
    axs[0][idx].set_title("psf_simga = %.1lf mm"%(psfsigma))
    axs[0][idx].set_xlabel("DETX (pixel)")
    axs[0][idx].set_ylabel("DETY (pixel)")
    
    
    xbins = np.linspace(np.floor(detx_crab+0.5) - delxysiz/2.0, np.floor(detx_crab+0.5) + delxysiz/2.0, delxysiz+1)
    delxhist, xedges = np.histogram(vdetx[vcut], xbins)
    axs[1][idx].stairs(delxhist, xedges, fill=True) #, color='skyblue', edgecolor='black')
    axs[1][idx].set_xlim(xedges[0], xedges[-1]) 


    axs[1][idx].set_xlabel("DETX (pixel)")
    axs[1][idx].set_ylabel("Counts")


plt.ion()
plt.show()
    
