#!/usr/bin/env python

import os, sys
import numpy as np
import astropy.io.fits as pyfits

import matplotlib.pyplot as plt
import matplotlib.colors as colors

### for detecor coordinate 
from teldef import Teldef
import telcoord


#fig, ax = plt.subplots(figsize=(5,5))
#ax.axis('equal')

fig2, ax2 = plt.subplots(figsize=(6,4))
ax2.set_yscale('log')
ax2.set_xscale('log')

#rspfname = "../refdata/hzx_erosita.rmf"
hzxsimdir = os.getenv("HZXSIMDIR")
rspfname = hzxsimdir+os.sep+"refdata/hzx_erosita.rmf"
hdul = pyfits.open(rspfname)
hdu = hdul["EBOUNDS"]
vemin = hdu.data.field("E_MIN")
vemax = hdu.data.field("E_MAX")
velh = np.append(vemin, vemax[-1:]) 
hdul.close()


### set up event data
obsid = "100000"  ### Crab
detid = 7

teldefverid = 2 ### teldef version
dirname = obsid 

phabinsize = 10
#phabinsize = 5
defbins = np.linspace(0.0, 1024.0, 1025)
newbins = defbins[::phabinsize]
velhbins = velh[::phabinsize]
vebinw   = velhbins[1:]-velhbins[:-1]

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


#### Crab cut
ra_crab  = 83.63321 
dec_crab = 22.01446
vcut_crab = np.logical_and(np.fabs(vra_obj-ra_crab)<0.0001, np.fabs(vra_obj-ra_crab)<0.0001)

### select for nxb, x-ray backround
vcut_nxb  = vphene<0
vcut_sxrb = np.logical_and(np.logical_not(vcut_nxb), np.logical_not(vcut_crab))  

## detector region cut
teldeffname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v%03d.fits"%(detid, teldefverid)
this_teldef = Teldef(teldeffname) 
detx_crab, dety_crab = telcoord.radec2detpix(ra_crab, dec_crab, attrot, this_teldef)
print("Crab position (DETX, DETY) = (%lf, %lf)"%(detx_crab, dety_crab))

#regrad = 10.0 ### pixel
#regrad = 20.0 ### pixel
regrad = 500.0 ### pixel, no cut
vcut_reg1 = (vdetx-detx_crab)**2 +(vdety-dety_crab)**2 < regrad**2  

vcut_crab = np.logical_and(vcut_crab, vcut_reg1)
vcut_nxb  = np.logical_and(vcut_nxb,  vcut_reg1)
vcut_sxrb = np.logical_and(vcut_sxrb, vcut_reg1)


### background
spechist, xedges = np.histogram(vpi[vcut_nxb], newbins)
vy  = spechist/vebinw/expos
vye = np.sqrt(spechist)/vebinw/expos
vx  = (velhbins[1:]+velhbins[:-1])/2.0
vxe = (velhbins[1:]-velhbins[:-1])/2.0
#ax2.stairs(vy, velhbins)
ax2.errorbar(vx, vy, xerr=vxe, yerr=vye, fmt='.', ms=3, lw=1) #velhbins)


### sxrb
spechist, xedges = np.histogram(vpi[vcut_sxrb], newbins)
#ax2.stairs(spechist/vebinw/expos, velhbins)
vy  = spechist/vebinw/expos
vye = np.sqrt(spechist)/vebinw/expos
#vx  = (velhbins[1:]+velhbins[:-1])/2.0
#vxe = (velhbins[1:]-velhbins[:-1])/2.0
#ax2.stairs(vy, velhbins)
ax2.errorbar(vx, vy, xerr=vxe, yerr=vye, fmt='.', ms=3, lw=1) #velhbins)


### crab
spechist, xedges = np.histogram(vpi[vcut_crab], newbins)
#ax2.stairs(spechist/vebinw/expos, velhbins)
vy  = spechist/vebinw/expos
vye = np.sqrt(spechist)/vebinw/expos
#vx  = (velhbins[1:]+velhbins[:-1])/2.0
#vxe = (velhbins[1:]-velhbins[:-1])/2.0
#ax2.stairs(vy, velhbins)
ax2.errorbar(vx, vy, xerr=vxe, yerr=vye, fmt='.', ms=3, lw=1, color='k') #velhbins)


ax2.set_xlim(0.1, 10)
ax2.set_xticks([0.1, 1, 10], ["0.1", "1", "10"])
ax2.set_xlabel("Energy (keV)")
ax2.set_ylabel(r"Counts keV$^{-1}$ s$^{-1}$")

plt.ion()
plt.show()

