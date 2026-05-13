#!/usr/bin/env python

import os, sys
import numpy as np

import astropy.io.fits as pyfits
from astropy.time import Time

import eptime

class ObsLog :
    def __init__(self, obslogfname) :


        hdulist = pyfits.open(obslogfname)
        hdu = hdulist[1]
        self.nrows = hdu.header["NAXIS2"]
        self.vobsid  = hdu.data.field("OBSID")
        self.vtstart = hdu.data.field("TSTART").astype('d')
        self.vtstop  = hdu.data.field("TSTOP").astype('d')
        self.vexpos  = hdu.data.field("EXPOSURE").astype('d')
        self.vqparam = hdu.data.field("QPARAMS")
        hdulist.close()

        return

    def get_obsid_param(self, obsid) :

        vidx = (self.vobsid == obsid)
        if np.any(vidx) : 
            tstart = self.vtstart[vidx]
            tstop  = self.vtstop[vidx]
            expos  = self.vexpos[vidx]
            qparam = self.vqparam[vidx]
        else :
            tstart = -1.0
            tstop  = -1.0
            expos  = -1.0
            qparam = []

        return tstart, tstop, expos, qparam 

    
    def get_paramlist(self, tstart, tstop) :
        vidx = np.all((tstart<=self.vtstart, self.vtstop<=tstop), axis=0)
        vobsid  = self.vobsid[vidx] 
        vtstart = self.vtstart[vidx] 
        vtstop  = self.vtstop[vidx]  
        vexpos  = self.vexpos[vidx] 
        vqparam = self.vqparam[vidx]
        return vobsid, vtstart, vtstop, vexpos, vqparam
        

if __name__=="__main__" :

    
    mjd = 60000
    mjdutc = Time(mjd, scale='utc', format='mjd')
    print(mjdutc)
    eptstart = eptime.mjdutc2eptime(mjd)
    print(eptstart)
    
    obslogfname = "../obsloglist.fits"

    obslog = ObsLog(obslogfname)
    
    #eptstart = ept
    eptstop  = eptstart+86400.
    vobsid, vtstart, vtstop, vexpos, vqparam = obslog.get_paramlist(eptstart, eptstop)
    nobs = len(vobsid)
    print("nobs = ", nobs)
    for idx in range(nobs) :
        print(vobsid[idx], vtstart[idx], vtstop[idx], vexpos[idx], vqparam[idx])
    
