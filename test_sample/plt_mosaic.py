#!/usr/bin/env python

import os, sys
import numpy as np
from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
import astropy.wcs as wcs
from astropy import units as u

import matplotlib.pyplot as plt
import matplotlib.colors as colors


import telcoord

### mosaic image dimension
sky_xsiz = 1200
sky_ysiz = 1200
xybinsize = 2

mosimg = np.zeros((sky_xsiz, sky_ysiz)) 

xbins = np.linspace(0.5, sky_xsiz+0.5, sky_xsiz+1)
ybins = np.linspace(0.5, sky_ysiz+0.5, sky_ysiz+1)


### TELDEF
#verid = 1
verid = 2
teldefdir = os.getenv("HZXSIMDIR")+os.sep+"refdata/teldef"
teldefver = '20251216v%03d'%(verid)


obsid = "100000"  ### Crab region
#obsid = "100001"  ### Galactic center
dirname = obsid 
#dirname = obsid+".def"

DRAW_FOV = True
#DRAW_FOV = False


### get wcs info
evtfname = dirname+os.sep+"hzx%sxm01.evt"%(obsid)
hdul = pyfits.open(evtfname)
hdu = hdul["EVENTS"]
q1 = hdu.header["QPARAM1"] 
q2 = hdu.header["QPARAM2"] 
q3 = hdu.header["QPARAM3"] 
q4 = hdu.header["QPARAM4"]

cdelt1 = hdu.header["TCDLT6"] ### DETX
cdelt2 = hdu.header["TCDLT7"] ### DETY
attquat = np.array([q1,q2,q3,q4])
print("attquat = ", attquat)
attrot = telcoord.attquat2rot(attquat)
radecroll = telcoord.attrot2radecroll(attrot, alignr=R.identity())
#radecroll = telcoord.skyrot2radecroll(attrot)
print(radecroll)
ra_cen  = radecroll[0]
dec_cen = radecroll[1]
        
wcs_mos = wcs.WCS(naxis=2)
wcs_mos.wcs.ctype = ["RA---TAN", "DEC--TAN"]
#wcs_mos.wcs.ctype = ["RA---CAR", "DEC--CAR"]
wcs_mos.wcs.crval = [ra_cen, dec_cen]
wcs_mos.wcs.cdelt = [cdelt1*xybinsize, cdelt2*xybinsize]
wcs_mos.wcs.crpix = [sky_xsiz/2., sky_ysiz/2.]
wcs_mos.wcs.pc = [[-1, 0], [0, 1]]        

hdul.close()



#fig, ax = plt.subplots(figsize=(6,6), subplot_kw=dict(projection=wcs_mos))
fig, ax = plt.subplots(figsize=(8,8), subplot_kw=dict(projection=wcs_mos))
#ax.axis('equal')
#ax.set_xlim(-50,1250)
#ax.set_ylim(-50,1250)

detid_list = range(1,17)
for detid in detid_list :
    #evtfname = "temp_xrbkg_xm%02d.evt"%(detid)
    evtfname = dirname+os.sep+"hzx%sxm%02d.evt"%(obsid, detid)
    hdul = pyfits.open(evtfname)
    hdu = hdul["EVENTS"]
        
    vra  = hdu.data.field("RA")
    vdec = hdu.data.field("DEC")
    vxpix, vypix = wcs_mos.wcs_world2pix(vra, vdec, 1)
    
    hist, xedges, yedges = np.histogram2d(vxpix, vypix, [xbins, ybins])
    mosimg += hist 
    #ax.plot(vxpix, vypix, '.', ms=1)
    ax.plot(vxpix, vypix, '.', ms=0.1, color='k')
    hdul.close()



ax.grid(color='k', ls='dotted')
ax.set(xlabel='Right Ascension', ylabel='Declination')

overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='magenta', ls='dotted', lw=2)
overlay[0].set_ticks([0.0, 180.]*u.degree)
overlay[1].set_ticks([0.0]*u.degree)
#overlay[0].set_axislabel('Galactic Longitude')
#overlay[1].set_axislabel('Galactic Latitude')


#figfname = dirname+"_raw.png"
#figfname = dirname+"a_raw.png"
#fig.savefig(figfname)

### FOV regions
from teldef import Teldef
from regions import Regions, PixelRegion, SkyRegion
from telfovtool import gen_telfovreg


if DRAW_FOV :
    teldef_list = []
    for detid in detid_list :
        teldeffname = teldefdir+os.sep+"hiz_xm%02d_teldef_%s.fits"%(detid, teldefver) 
        teldef = Teldef(teldeffname)
        teldef_list.append(teldef)


    ### create FOV region file
    regfname = "temp_fovskyimg.reg"
    color_list = ["black", "red", "forestgreen", "blue"]
    coordsys = 'icrs'
    gen_telfovreg(regfname, attrot, teldef_list, coordsys, color_list)
    
    regions = Regions.read(regfname, format='ds9')
    for i, region in enumerate(regions):
        print(i, region)
        pixreg = None
        if isinstance(region, SkyRegion) :
            pixreg = region.to_pixel(wcs_mos)
        elif isinstance(region, PixelRegion) :
            pixreg = region
        if pixreg is not None :
            pixreg.plot(ax=ax)
            #pixreg.plot(ax=ax2)
    
    #figfname = dirname+"a_fov.png"
    #fig.savefig(figfname)

plt.ion()
plt.show()

