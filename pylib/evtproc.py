#!/usr/bin/env python

import os
import numpy as np
from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits
import astropy.wcs as wcs
from astropy.time import Time

#import healpy


from teldef import Teldef
import telcoord
from mstime import MSTime, mstime_start, mstime_stop


class EventProc :
    def __init__(self) :
        ### default 
        self.teldef = None
        return
        
    #def setObsParams(self, teldef, tstart, tstop, qparam, objectName="XRAYSKY", obsid=None, inv=True) :
    ### Change default inv to False
    def setObsParams(self, teldef, tstart, tstop, qparam, objectName="XRAYSKY", obsid=None, inv=False) :

        self.teldef = teldef
        self.tstart = tstart
        self.tstop  = tstop
        self.qparam = qparam
        self.obsid  = obsid
        if inv :
            self.attrot = R.from_quat(qparam).inv()
        else :
            self.attrot = R.from_quat(qparam)

        self.object = objectName

        #### skyref
        self.skyref = telcoord.attrot2radecroll(self.attrot, self.teldef)
        
        
        return 


    def setStandardHeaderKeywords(self, hdu, extname) :

        telapse = self.tstop - self.tstart
        exposure = self.tstop - self.tstart

        utc_start = Time(MSTime.mstime2mjdutc(self.tstart), scale='utc', format='mjd')
        utc_stop  = Time(MSTime.mstime2mjdutc(self.tstop), scale='utc', format='mjd')

        
        hdu.header["EXTNAME" ] = (extname, "Extenstion name")
        hdu.header["TELESCOP"] = (self.teldef.telescop, "Telescope name")
        hdu.header["INSTRUME"] = (self.teldef.instrume, "Instrument name")
        hdu.header["DETNAM"  ] = (self.teldef.detnam  , "Detector Name")
        hdu.header["DATAMODE"] = ("PHOTON"       , "Instrument operating mode")
        hdu.header["TIMESYS" ] = ("TT"           , "Time system")
        hdu.header["MJDREFI" ] = (MSTime.mjdrefi, "Reference MJD Integer part")
        hdu.header["MJDREFF" ] = (MSTime.mjdreff, "Reference MJD fractional")
        hdu.header["TIMEREF" ] = ("LOCAL"       , "Reference time")
        hdu.header["TASSIGN" ] = ("SATELLITE"   , "Time assigned by clock")
        hdu.header["TIMEUNIT"] = ("s", "Time unit")
        hdu.header["TSTART"  ] = (self.tstart, "Start time")
        hdu.header["TSTOP"   ] = (self.tstop , "Stop time")
        hdu.header["TELAPSE" ] = (telapse    , "Elapsed time")
        hdu.header["EXPOSURE"] = (exposure   , "Exposure time")

        hdu.header["DATE-OBS"] = utc_start.datetime.strftime("%Y-%m-%d")
        hdu.header["TIME-OBS"] = utc_start.datetime.strftime("%H:%M:%S")
        hdu.header["DATE-END"] = utc_stop.datetime.strftime("%Y-%m-%d")
        hdu.header["TIME-END"] = utc_stop.datetime.strftime("%H:%M:%S")

        hdu.header["RADECSYS"] = ("FK5", "World Coordinate System")
        hdu.header["EQUINOX"]  = (2000., "Equinox for coordinate system")
        

        ### OBSID
        if self.obsid is not None :
            hdu.header["OBS_ID"] = self.obsid
        
        
        ### OBJECT Name
        hdu.header["OBJECT"] = (self.object, "Name of object")

        ### QPARAM
        hdu.header["QPARAM1"] = self.qparam[0]
        hdu.header["QPARAM2"] = self.qparam[1]
        hdu.header["QPARAM3"] = self.qparam[2]
        hdu.header["QPARAM4"] = self.qparam[3]

        ### SKYREF parameter
        hdu.header["RA_PNT"]  = self.skyref[0]
        hdu.header["DEC_PNT"] = self.skyref[1]
        hdu.header["PA_PNT"]  = self.skyref[2]


        
    def setEventHeaderKeywords(self, hdu, wcs_raw, wcs_det, wcs_sky) :

        teldef = self.teldef
        ### XY column
        hdu.header["TLMIN2"] = teldef.skyxpix1
        hdu.header["TLMAX2"] = teldef.sky_xsiz
        hdu.header["TLMIN3"] = teldef.skyypix1
        hdu.header["TLMAX3"] = teldef.sky_ysiz
    
        hdu.header["TCTYP2"] = wcs_sky.wcs.ctype[0]
        hdu.header["TCTYP3"] = wcs_sky.wcs.ctype[1]
        #hdu.header["TCDLT2"] = wcs_sky.wcs.cdelt[0]
        hdu.header["TCDLT2"] = -wcs_sky.wcs.cdelt[0] ### modified for xselect
        hdu.header["TCDLT3"] = wcs_sky.wcs.cdelt[1]
        hdu.header["TCRPX2"] = wcs_sky.wcs.crpix[0]
        hdu.header["TCRPX3"] = wcs_sky.wcs.crpix[1]
        hdu.header["TCRVL2"] = wcs_sky.wcs.crval[0]
        hdu.header["TCRVL3"] = wcs_sky.wcs.crval[1]

        ### modified below for xselect
        #hdu.header["TP2_2"]  = wcs_sky.wcs.pc[0][0] 
        #hdu.header["TP2_3"]  = wcs_sky.wcs.pc[0][1]
        #hdu.header["TP3_2"]  = wcs_sky.wcs.pc[1][0] 
        #hdu.header["TP3_3"]  = wcs_sky.wcs.pc[1][1] 


        ### RAWXY column
        hdu.header["TLMIN4"] = teldef.rawxpix1     
        hdu.header["TLMAX4"] = teldef.raw_xsiz     
        hdu.header["TLMIN5"] = teldef.rawypix1     
        hdu.header["TLMAX5"] = teldef.raw_ysiz     
                                                      
        hdu.header["TCTYP4"] = wcs_raw.wcs.ctype[0]
        hdu.header["TCTYP5"] = wcs_raw.wcs.ctype[1]
        hdu.header["TCDLT4"] = wcs_raw.wcs.cdelt[0]
        hdu.header["TCDLT5"] = wcs_raw.wcs.cdelt[1]
        hdu.header["TCRPX4"] = wcs_raw.wcs.crpix[0]  
        hdu.header["TCRPX5"] = wcs_raw.wcs.crpix[1]  
        hdu.header["TCRVL4"] = wcs_raw.wcs.crval[0]
        hdu.header["TCRVL5"] = wcs_raw.wcs.crval[1]
        
        hdu.header["TP4_4"]  = wcs_raw.wcs.pc[0][0]
        hdu.header["TP4_5"]  = wcs_raw.wcs.pc[0][1]
        hdu.header["TP5_4"]  = wcs_raw.wcs.pc[1][0]
        hdu.header["TP5_5"]  = wcs_raw.wcs.pc[1][1]


        ### DETXY column
        ### temporally copy of RAWXY image
        hdu.header["TLMIN6"] = teldef.rawxpix1     
        hdu.header["TLMAX6"] = teldef.raw_xsiz     
        hdu.header["TLMIN7"] = teldef.rawypix1     
        hdu.header["TLMAX7"] = teldef.raw_ysiz     
        
        hdu.header["TCTYP6"] = wcs_det.wcs.ctype[0]
        hdu.header["TCTYP7"] = wcs_det.wcs.ctype[1]
        hdu.header["TCDLT6"] = wcs_det.wcs.cdelt[0]
        hdu.header["TCDLT7"] = wcs_det.wcs.cdelt[1]
        hdu.header["TCRPX6"] = wcs_det.wcs.crpix[0]  
        hdu.header["TCRPX7"] = wcs_det.wcs.crpix[1]  
        hdu.header["TCRVL6"] = wcs_det.wcs.crval[0]
        hdu.header["TCRVL7"] = wcs_det.wcs.crval[1]
        
        hdu.header["TP6_6"]  = wcs_det.wcs.pc[0][0]
        hdu.header["TP6_7"]  = wcs_det.wcs.pc[0][1]
        hdu.header["TP7_6"]  = wcs_det.wcs.pc[1][0]
        hdu.header["TP7_7"]  = wcs_det.wcs.pc[1][1]


        ### PHA/PI column
        hdu.header["TLMIN10"] = 0     ### PHA
        hdu.header["TLMAX10"] = 1023
        hdu.header["TLMIN11"] = 0     ### PI
        hdu.header["TLMAX11"] = 1023


        
    def procEventToFile(self, ra, dec, vtime, vdetxpix, vdetypix, vpha, vene, outfname) :
    
        print("Processing simulated events")
    
        nevts = len(vtime)

        teldef = self.teldef
        attrot = self.attrot
        skyref = self.skyref
        
        #teldef = self.teldef
        #ra_opt, dec_opt, rotrad = telcoord.get_wcsprm_detalign(attrot, teldef.alignrot)
        #ra_opt, dec_opt, rotrad = telcoord.get_wcsprm_detalign(attrot, teldef.alignrot, vdet=np.array([0,0,1]).astype('d'))
        #skyref = telcoord.get_skyref(attrot, teldef)
        ra_opt  = skyref[0]
        dec_opt = skyref[1]
        rot_opt = skyref[2]
        rotrad  = np.radians(rot_opt)
        
        cos_r = np.cos(rotrad)
        sin_r = np.sin(rotrad)
        #print("RA, Dec, Rot =", ra_opt, dec_opt, np.degrees(rotrad))
        print("RA, Dec, Rot =", ra_opt, dec_opt, rot_opt)

        if isinstance(ra, np.ndarray) and isinstance(dec, np.ndarray) :
            vra_obj  = ra
            vdec_obj = dec 
        else :
            vra_obj  = np.zeros(nevts)+ra
            vdec_obj = np.zeros(nevts)+dec
            
    
        ### initialize WCS object
        wcs_raw = wcs.WCS(naxis=2)
        wcs_det = wcs.WCS(naxis=2)
        wcs_sky = wcs.WCS(naxis=2)
        
        cdeltx = np.degrees(np.arctan(teldef.det_xscl / teldef.focallen))
        cdelty = np.degrees(np.arctan(teldef.det_yscl / teldef.focallen))
    
        ### wcs for raw, det coordiantes
        for w in [wcs_raw, wcs_det, wcs_sky] :
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            w.wcs.crval = [ra_opt, dec_opt]
            w.wcs.cdelt = [cdeltx, cdelty]
            w.wcs.crpix = [teldef.optaxisx, teldef.optaxisy]
            w.wcs.pc = [[-cos_r, sin_r], [sin_r, cos_r]]
    
        ### wcs for sky coordiantes
        wcs_sky.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        wcs_sky.wcs.crval = [ra_opt, dec_opt]
        wcs_sky.wcs.cdelt = [cdeltx, cdelty]
        wcs_sky.wcs.crpix = [teldef.sky_xsiz/2., teldef.sky_ysiz/2.]
        wcs_sky.wcs.pc = [[-1, 0], [0, 1]]
        
        
        ### Event proces
        vdet = telcoord.detpix2vdet(vdetxpix, vdetypix, teldef)
        #print("vdet = ", vdet)
        vra, vdec = telcoord.vdet2radec(vdet, attrot, teldef.alignrot)
        #print("vra, vdec = ", vra, vdec)
        vxpix, vypix = wcs_sky.wcs_world2pix(vra, vdec, 1)
    
        ### Event HDU
        ### sort array with time
        vsort = np.argsort(vtime)
        
        col1 = pyfits.Column(name="TIME", format="D", unit="s", array=vtime[vsort])
    
        col2 = pyfits.Column(name="X",    format="I", unit="pixel", array=vxpix[vsort]) 
        col3 = pyfits.Column(name="Y",    format="I", unit="pixel", array=vypix[vsort]) 
    
        #col4 = pyfits.Column(name="RAWX", format="I", unit="pixel", array=vrawx)
        #col5 = pyfits.Column(name="RAWY", format="I", unit="pixel", array=vrawy)
        col4 = pyfits.Column(name="RAWX", format="I", unit="pixel", array=vdetxpix[vsort])
        col5 = pyfits.Column(name="RAWY", format="I", unit="pixel", array=vdetypix[vsort])
    
        col6 = pyfits.Column(name="DETX", format="I", unit="pixel", array=vdetxpix[vsort]) ### copy raw to det
        col7 = pyfits.Column(name="DETY", format="I", unit="pixel", array=vdetypix[vsort])
    
        col8 = pyfits.Column(name="RA",   format="E", unit="degree", array=vra[vsort])
        col9 = pyfits.Column(name="DEC",  format="E", unit="degree", array=vdec[vsort])
    
        col10 = pyfits.Column(name="PHA", format="I", unit="channel", array=vpha[vsort])  ### copy pha to pi
        col11 = pyfits.Column(name="PI",  format="I", unit="channel", array=vpha[vsort])
    
        #col12 = pyfits.Column(name="RA_OBJ",  format="E", unit="degree", array=np.ones(nevts,'f')*ra_obj)
        #col13 = pyfits.Column(name="DEC_OBJ", format="E", unit="degree", array=np.ones(nevts,'f')*dec_obj)
        col12 = pyfits.Column(name="RA_OBJ",  format="E", unit="degree", array=vra_obj[vsort])
        col13 = pyfits.Column(name="DEC_OBJ", format="E", unit="degree", array=vdec_obj[vsort])
        col14 = pyfits.Column(name="PHOTON_E", format="E", unit="keV", array=vene[vsort])
        #col14 = pyfits.Column(name="PHOTON_E", format="D", unit="keV", array=vene[vsort])
        

        cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14])
        evthdu = pyfits.BinTableHDU.from_columns(cols)
        
        extname = "EVENTS"
        self.setStandardHeaderKeywords(evthdu, extname)
        self.setEventHeaderKeywords(evthdu, wcs_raw, wcs_det, wcs_sky)

        ### GTI HDU
        col_start = pyfits.Column(name="START", format="D", unit="s", array=np.array([self.tstart]))
        col_stop  = pyfits.Column(name="STOP",  format="D", unit="s", array=np.array([self.tstop]))
        cols   = pyfits.ColDefs([col_start, col_stop])
        gtihdu = pyfits.BinTableHDU.from_columns(cols)
        
        extname = "GTI"
        self.setStandardHeaderKeywords(gtihdu, extname)

        ### output file
        prihdr = pyfits.Header() 
        prihdu = pyfits.PrimaryHDU(header=prihdr)
        
        fptr_out = pyfits.HDUList([prihdu, evthdu, gtihdu])
        fptr_out.writeto(outfname)
    
        return
 
