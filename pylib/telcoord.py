#!/usr/bin/env python

import os
import math
import numpy as np
from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
import astropy.wcs as wcs

from teldef import Teldef


def rawpix2detpix(rawxpix, rawypix) :
    ### correction of optical distortion
    ### temporally copy rawxy to detxy
    detxpix = rawxpix
    detypix = rawypix
    return detxpix, detypix


def detpix2detmm(detxpix, detypix, teldef) :
    detxmm =  (detxpix - teldef.detxpix_cen) * teldef.det_xscl 
    detymm =  (detypix - teldef.detypix_cen) * teldef.det_yscl
    return detxmm, detymm 
    
def detmm2detpix(detxmm, detymm, teldef) :
    detxpix =  detxmm / teldef.det_xscl + teldef.detxpix_cen
    detypix =  detymm / teldef.det_yscl + teldef.detypix_cen
    return detxpix, detypix


def detmm2vdet(detxmm, detymm, teldef) :
    #if isinstance(detxmm, np.ndarray) :
    #    nelem = len(detxmm)  ### get nelements
    #    vdetx = np.zeros(nelem).astype('d') + teldef.focallen
    #else :
    #    vdetx = teldef.focallen
    #vdety = -detxmm
    #vdetz =  detymm
    
    if isinstance(detxmm, np.ndarray) :
        nelem = len(detxmm)  ### get nelements
        vdetz = np.zeros(nelem).astype('d') + teldef.focallen
    else :
        vdetz = +teldef.focallen
    vdetx = -detxmm  ### detx sign is swapped. 
    vdety =  detymm
    
    vdet = np.array([vdetx, vdety, vdetz]) 
    return vdet


def vdet2detmm(vdet, teldef) :

    #vdety = vdet[1] * teldef.focallen/vdet[0] 
    #vdetz = vdet[2] * teldef.focallen/vdet[0] 
    #detxmm = -vdety
    #detymm =  vdetz

    vdetx = vdet[0] * teldef.focallen/vdet[2] 
    vdety = vdet[1] * teldef.focallen/vdet[2] 
    detxmm = -vdetx  ### detx sign is swapped. 
    detymm =  vdety


    return detxmm, detymm
    


def detpix2vdet(detxpix, detypix, teldef) :

    detxmm, detymm = detpix2detmm(detxpix, detypix, teldef)
    vdet = detmm2vdet(detxmm, detymm, teldef)
    return vdet
    


def vdet2detpix(vdet, teldef) :

    detxmm, detymm = vdet2detmm(vdet, teldef)
    detxpix, detypix = detmm2detpix(detxmm, detymm, teldef)
    return detxpix, detypix
    


def vdet2vsky(vdet, attrot, alignrot, skyrot=None) :

    if len(vdet.shape)==1 :
        naxis2 = 0
        if skyrot is None :
            vsat = alignrot.apply(vdet, inverse=True) ### CALDB teldef matrix 
            #vsky = attrot.apply(vsat, inverse=True)  ### Old attrot definition
            vsky = attrot.apply(vsat)
        else :
            #vsky = skyrot.apply(vdet, inverse=True)  ### Old attrot definition
            vsky = skyrot.apply(vdet)
    else :
        naxis2 = vdet.shape[1]
        vvdet = vdet.transpose()
        if skyrot is None :
            vvsat = alignrot.apply(vvdet, inverse=True) ### CALDB teldef matrix 
            #vvsky = attrot.apply(vvsat, inverse=True)  ### Old attrot definition
            vvsky = attrot.apply(vvsat)  
        else :
            #vvsky = skyrot.apply(vvdet, inverse=True)   ### Old attrot definition
            vvsky = skyrot.apply(vvdet)  

        vsky = vvsky.transpose()
    return vsky



def vsky2vdet(vsky, attrot, alignrot, skyrot=None) :
    if len(vsky.shape)==1 :
        if skyrot is None :
            #vsat  = attrot.apply(vsky)  ### Old attrot definition
            vsat  = attrot.apply(vsky, inverse=True)
            vdet  = alignrot.apply(vsat)  ### CALDB teldef matrix
        else:
            vdet = skyrot.apply(vsky)
        
    else :
        vvsky = vsky.transpose()
        if skyrot is None :
            #vvsat = attrot.apply(vvsky) ### Old attrot definition
            vvsat = attrot.apply(vvsky, inverse=True) 
            vvdet = alignrot.apply(vvsat)  ### CALDB teldef matrix
        else :
            vvdet = skyrot.apply(vvsky)
            
        vdet = vvdet.transpose()
    return vdet


def detpix2vsky(detxpix, detypix, attrot, teldef) :
    vdet = detpix2vdet(detxpix, detypix, teldef)
    vsky = vdet2vsky(vdet, attrot, teldef.alignrot) 
    return vsky



def vsky2detpix(vsky, attrot, teldef) :
    vdet = vsky2vdet(vsky, attrot, teldef.alignrot)
    detxpix, detypix = vdet2detpix(vdet, teldef)
    return detxpix, detypix



def radec2detpix(ra, dec, attrot, teldef) :
    vdet = radec2vdet(ra, dec, attrot, teldef.alignrot)
    detxpix, detypix = vdet2detpix(vdet, teldef)
    return detxpix, detypix


def detpix2radec(detxpix, detypix, attrot, teldef) :
    vdet = detpix2vdet(detxpix, detypix, teldef)
    ra, dec = vdet2radec(vdet, attrot, teldef.alignrot)
    return ra, dec



def vsky2radec(vsky) :
    #print(vsky, vsky.shape, len(vsky.shape))
    #if vsky.shape
    vskyx = vsky[0]
    vskyy = vsky[1]
    vskyz = vsky[2]
    
    r    = np.sqrt(vskyx**2+vskyy**2+vskyz**2)
    dec  = np.degrees(np.arcsin(vskyz/r))       
    ra   = np.degrees(np.arctan2(vskyy, vskyx)) 

    ### add 360. to negative values
    ra += (ra<0.)*360.
    
    return ra, dec


def radec2vsky(ra, dec) :

    vskyx = np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    vskyy = np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    vskyz = np.sin(np.radians(dec))

    vsky  = np.array((vskyx, vskyy, vskyz))
    
    return vsky


def vdet2radec(vdet, attrot, alignrot) :
    vsky = vdet2vsky(vdet, attrot, alignrot)
    ra, dec = vsky2radec(vsky)
    return ra, dec


def radec2vdet(ra, dec, attrot, alignrot) :
    vsky = radec2vsky(ra, dec) 
    vdet = vsky2vdet(vsky, attrot, alignrot)
    return vdet



def euler2radecroll(euler, roff=0.0, rsign=-1.0) :

    ra   = euler[0]
    dec  = 90.-euler[1]
    roll = euler[2]-90.

    if ra<0.0 : ra+=360.0
    
    roll = (roll + roff) * rsign
    if roll<0.0 : roll+=360.0

    skyref = np.array((ra, dec, roll))
    return skyref



def radecroll2euler(skyref, roff=0.0, rsign=-1.0) :

    euler = np.zeros(3)
    euler[0] = skyref[0]     # ra to euler1
    euler[1] = 90.-skyref[1] # dec to euler2
    euler[2] = skyref[2]*rsign + 90.- roff ### roll to euler3

    return euler



#def skyref2rot(skyref, roff=0.0, rsign=-1.0, as_euler=False) :
def radecroll2skyrot(skyref, roff=0.0, rsign=-1.0, as_euler=False) :

    if as_euler :
        euler = skyref
    else :
        #euler = np.zeros(3)
        #euler[0] = skyref[0]     # ra to euler1
        #euler[1] = 90.-skyref[1] # dec to euler2
        #euler[2] = skyref[2]*rsign + 90.- roff ### roll to euler3
        euler = radecroll2euler(skyref, roff, rsign)
        
    skyrot = atteuler2rot(euler)
    return skyrot



### calculate ra, dec, roll from attitude and teldef alignment parameters 
### codes are based on heasoft attitude libraries
### in attitude/lib/coordfits/align.c
def skyrot2radecroll(skyrot, roff=0.0, rsign=-1.0, as_euler=False) :

    ### !!! Bug !!!
    ### skyrot = alignrot * attrot             : vsky -> vdet
    ### skyrot.inv = attrot.inv * alingot.inv  : vdet -> vsky

    ### Correct
    ### skyrot = attrot * alignrot.inv()       : vsky -> vdet
    ### skyrot.inv = alingot * attrot.inv()    : vdet -> vsky

    #euler = skyrot.inv().as_euler('ZYZ', degrees=True) ### FTOOLS attitude definition
    euler = attrot2euler(skyrot) 

    if as_euler :
        return euler

    #ra  = euler[0]
    #dec  = 90.-euler[1]
    #roll = euler[2]-90.
    #
    #if ra<0.0 : ra+=360.0
    #
    #roll = (roll + roff) * rsign
    #if roll<0.0 : roll+=360.0
    #
    #skyref = np.array((ra, dec, roll))

    skyref = euler2radecroll(euler, roff, rsign)
    
    return skyref



### calculate ra, dec, roll from attitude and teldef alignment parameters 
### codes are based on heasoft attitude libraries
### in attitude/lib/coordfits/align.c
def attrot2radecroll(attrot, teldef=None,
                     alignr=R.from_matrix([[0,1,0],[0,0,1],[1,0,0]]), roff=0.0, rsign=-1.0,
                     as_euler=False) :
    
    if teldef == None :
        ### default
        alignrot = alignr
        rolloff  = roff
        rollsign = rsign
    else :
        alignrot = teldef.alignrot
        rolloff  = teldef.rolloff
        rollsign = teldef.rollsign

    #skyrot =  alignrot * attrot  ### original, bug?
    ### test
    skyrot =  attrot * alignrot.inv() ### looks good

    return skyrot2radecroll(skyrot, rolloff, rollsign, as_euler) 
    


def radecroll2attrot(skyref, teldef=None,
                     alignr=R.from_matrix([[0,1,0],[0,0,1],[1,0,0]]), roff=0.0, rsign=-1.0,
                     as_euler=False) :
    
    if teldef == None :
        ### default
        alignrot = alignr
        rolloff  = roff
        rollsign = rsign
    else :
        alignrot = teldef.alignrot
        rolloff  = teldef.rolloff
        rollsign = teldef.rollsign

    
    skyrot = radecroll2skyrot(skyref, rolloff, rollsign, as_euler)

    ### attrot = alignrot.inv() * skyrot ### original, bug
    attrot = skyrot * alignrot  ### should be ok

    return attrot



def attquat2rot(quat) :
    ### Definition of attitude quaternion to rotation function 
    #attrot = R.from_quat(quat).inv()  ### Old attrot definition
    attrot = R.from_quat(quat) 
    return attrot


def attrot2quat(attrot) :
    #quat = attrot.inv().as_quat()  ### Old attrot definition 
    quat = attrot.as_quat()
    return quat


def atteuler2rot(euler) :
    #attrot = R.from_euler('ZYZ', euler, degrees=True).inv()  ### Old attrot definition
    attrot = R.from_euler('ZYZ', euler, degrees=True)  
    return attrot


def attrot2euler(attrot) :
    #euler = attrot.inv().as_euler('ZYZ', degrees=True)  ### Old attrot definition
    euler = attrot.as_euler('ZYZ', degrees=True)
    return euler



def calcXYOffAxis(ra, dec, attrot, alignrot) :
    vdet = radec2vdet(ra, dec, attrot, alignrot)
    xx = -vdet[0]
    yy =  vdet[1]
    zz =  vdet[2]
    theta_x = np.degrees(np.arccos(zz/np.sqrt(xx*xx+zz*zz)))
    theta_y = np.degrees(np.arccos(zz/np.sqrt(yy*yy+zz*zz)))
    return theta_x, theta_y



if __name__=="__main__" :


    #qparam = [0,0,0,1]
    #qparam = [-0.70955226, -0.00353587,  0.00138493,  0.70464259] ### bad
    qparam = [-0.58883817,  0.00332964, -0.00195868,  0.80824173]
    #attrot = R.from_quat(qparam).inv()
    attrot = attquat2rot(qparam)
    print(attrot.as_matrix())

    #detid = 1

    detid = 7
    #detid = 8

    #detid = 22
    
    ### CALDB teldef
    #caldbdir = os.getenv("CALDB")
    #teldeffname = "hizx_teldef_temp.fits"
    teldeffname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v001.fits"%(detid)
    teldef = Teldef(teldeffname)
    #print(teldef.alignmat)

    ### test vdet2radec
    vdet = np.array([[0], [0], [1]]).astype('d')
    #print(np.matmul(teldef.alignmat, vdet))

    vra, vdec = vdet2radec(vdet, attrot, teldef.alignrot)
    print("RA, Dec = ", vra[0], vdec[0])

    #ra, dec, rotrad =  get_wcsprm_detalign(quat, teldef.rm_det2sat)

    ra_opt, dec_opt, rotdeg =  attrot2radecroll(attrot, teldef)
    print("ra_pot, dec_opt, rotdeg = ", ra_opt, dec_opt, rotdeg)
    #ra_opt  = vra_opt[0]
    #dec_opt = vdec_opt[0]
    #rotrad  = vrotrad[0]
    cos_r = np.cos(np.radians(rotdeg))
    sin_r = np.sin(np.radians(rotdeg))
    print (ra_opt, dec_opt, cos_r, sin_r)


    detx0 = teldef.detxpix1  ### 1st pixel
    detx1 = teldef.detxpix1 + teldef.det_xsiz  ### last pixel
    dety0 = teldef.detypix1
    dety1 = teldef.detypix1 + teldef.det_ysiz

    vxpix = np.array([detx0, detx1, detx0, detx1]).astype('d')
    vypix = np.array([dety0, dety0, dety1, dety1]).astype('d')
    num = len(vxpix)


    print("\n### Check attitude, teldef functions ###")
    vra, vdec = detpix2radec(vxpix, vypix, attrot, teldef)
    vxpix2, vypix2 = radec2detpix(vra, vdec, attrot, teldef)
    for idx in range(num) :
        print("%6.1lf, %6.1lf ->  %14.10lf, %14.10lf -> %7.2lf, %7.2lf  "%(vxpix[idx], vypix[idx], vra[idx], vdec[idx], vxpix2[idx], vypix2[idx]))


    ### create wcs
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [teldef.detxpix_cen, teldef.detypix_cen]
    #w.wcs.crval = [ra_opt, dec_opt]
    w.wcs.crval = [ra_opt, dec_opt]
    cdeltx = np.degrees(np.arctan(teldef.det_xscl / teldef.focallen))
    cdelty = np.degrees(np.arctan(teldef.det_yscl / teldef.focallen))
    w.wcs.cdelt = [cdeltx, cdelty]

    w.wcs.pc = [[-cos_r, sin_r], [sin_r, cos_r]]
    ###
    ###  same as the Haiwu's code.
    ###  [[-1 0]  * [[cos_r -sin_r]  = [[-cos_r sin_r] 
    ###   [ 0 1]]    [sin_r  cos_r]]     [sin_r cos_r]]
    ###
    
    print(w)

    print("\n### Check wcs functions ###")
    vra, vdec = w.wcs_pix2world(vxpix, vypix, 1)
    vxpix2, vypix2 = w.wcs_world2pix(vra, vdec, 1)
    for idx in range(num) :
        print("%6.1lf, %6.1lf ->  %14.10lf, %14.10lf -> %7.2lf, %7.2lf"%(vxpix[idx], vypix[idx], vra[idx], vdec[idx], vxpix2[idx], vypix2[idx]))



    #### Check Crab position
    print("\n### Check Crab position ###")
    ra_obj  = 83.633
    dec_obj = 22.014

    print("RA, Dec = ", ra_obj, dec_obj)
    #print("### radec2detpix ###")
    detxpix, detypix = radec2detpix(ra_obj, dec_obj, attrot, teldef)
    print("radec2detpix: DETXPIX, DETYPIX = %6.2lf, %6.2lf"%(detxpix, detypix))
    
    detxpix2, detypix2 = w.wcs_world2pix(ra_obj, dec_obj, 1)
    print("wcs_world2pix: DETXPIX, DETYPIX = %6.2lf, %6.2lf"%(detxpix2, detypix2))
        
    
