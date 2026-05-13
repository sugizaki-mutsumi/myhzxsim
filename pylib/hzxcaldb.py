#!/usr/bin/env python

import os
import numpy as np
#from scipy.spatial.transform import Rotation as R

import astropy.io.fits as pyfits


#import datetime
from astropy.time import Time

#import eptime

from teldef import Teldef

from hzxpsf import HZXPSF
from hzxvign import HZXVign


from respmatrix import ResponseMatrix, AncillaryResponse



def get_caldb_fname(caldbidxfname, telescop, instrume, detnam, filterstr, codenam, fnamesub="") :

    fnamelist = []
    if not os.path.isfile(caldbidxfname) :
        return fnamelist

    hdul = pyfits.open(caldbidxfname)
    hdu = hdul["CIF"]
    nrows = hdu.header["NAXIS2"]

    idxtable = hdu.data
    vtelescop = idxtable.field("TELESCOP")
    vinstrume = idxtable.field("INSTRUME")
    vdetnam   = idxtable.field("DETNAM")
    vfilter   = idxtable.field("FILTER")
    vcal_dir  = idxtable.field("CAL_DIR")
    vcal_file = idxtable.field("CAL_FILE")
    vcal_cnam = idxtable.field("CAL_CNAM")
    vcal_qual = idxtable.field("CAL_QUAL")
    
    #vgood = np.all([(vtelescop==telescop), (vinstrume==instrume), (vdetnam==detnam), (vcal_cnam==codenam), np.char.find(vcal_file, fnamesub)>=0], axis=0)  

    vgood1 = np.all([vtelescop==telescop, vinstrume==instrume, vdetnam==detnam, vcal_cnam==codenam, np.char.find(vcal_file, fnamesub)>=0], axis=0)  
    vgood2 = np.all([vgood1, vcal_qual==0], axis=0)
    ngood1 = (vgood1*1).sum()
    ngood2 = (vgood2*1).sum()

    if ngood2>0 :
        vgood = vgood2
    elif ngood1>0 :
        vgood = vgood1
    else :
        vgood == None
        
    if vgood is not None :
        fnamelist = np.char.add(np.char.strip(vcal_dir[vgood]), os.sep)
        fnamelist = np.char.add(fnamelist, np.char.strip(vcal_file[vgood]))

    return fnamelist




class HzxCaldb :
    #def __init__(self, cmosid=1, teldefname="CALDB", rmffname="CALDB", arffname="CALDB", vigfname="CALDB", caldbdir="CALDB") :
    def __init__(self, detid=1, teldefname="CALDB", rmffname="CALDB", arffname="CALDB", vigfname="CALDB", caldbdir="CALDB") :
        
        ### CALDB mission/instrument name
        mission   = "HIZ"  ### CALDB mission directory name
        subsystem = "EAGLE" ### CALDB subsystem directory name
        telescop  = "HIZ"
        instrume  = "EAGLE"
        filterstr = "-" ### dummy

        
        self.detid = detid
        detnam = "XM%02d"%(detid)
        
        if caldbdir=="CALDB" :
            caldbdir = os.getenv("CALDB") 
        caldbinstdir = caldbdir+os.sep+"data/%s/%s"%(mission, subsystem)
        caldbidxfname = caldbinstdir+os.sep+"caldb.indx"

        
        ####  TELDEF 
        if teldefname=="CALDB" :
            codenam  = "TELDEF" 
            calfname_list = get_caldb_fname(caldbidxfname, telescop, instrume, detnam, filterstr, codenam)
            ncalfiles = len(calfname_list) 
            if len(calfname_list)==0 :
                self.teldeffname = None
            else :
                #self.teldeffname = caldbdir+os.sep+calfname_list[0]
                self.teldeffname = caldbdir+os.sep+calfname_list[-1]
        elif teldefname=="NONE" :
            self.teldeffname = None
        else :
            self.teldeffname = teldefname

        if self.teldeffname != None :
            print("Teldef file = %s"%(self.teldeffname))
            self.teldef = Teldef(self.teldeffname)

            
        ### Response Matrix
        if rmffname=="CALDB" :
            codenam  = "MATRIX" 
            calfname_list = get_caldb_fname(caldbidxfname, telescop, instrume, detnam, filterstr, codenam, fnamesub)
            ncalfiles = len(calfname_list) 
            if len(calfname_list)==0 :
                self.rmffname =  None
            else :
                #self.rmffname = caldbdir+os.sep+calfname_list[0]
                self.rmffname = caldbdir+os.sep+calfname_list[-1]
        elif rmffname=="NONE" :
            self.rmffname =  None
        else :
            self.rmffname = rmffname

        if self.rmffname!=None :
            print("RMF file = %s"%(self.rmffname))
            self.rmf = ResponseMatrix(self.rmffname)


        ### Ancillary reponnse file
        if arffname=="CALDB" :
            codenam  = "SPECRESP" 
            calfname_list = get_caldb_fname(caldbidxfname, telescop, instrume, detnam, filterstr, codenam, fnamesub)
            ncalfiles = len(calfname_list) 
            if len(calfname_list)==0 :
                self.arffname = None
            else :
                #self.arffname = caldbdir+os.sep+calfname_list[0]
                self.arffname = caldbdir+os.sep+calfname_list[-1]
        elif arffname=="NONE" :
            self.arffname = None
        else :
            self.arffname = arffname
            
        if self.arffname!=None :
            print("ARF file = %s"%(self.arffname))
            self.arf = AncillaryResponse(self.arffname)
            

        ### Vignetting
        if vigfname=="CALDB" :
            codenam  = "VIGNET" 
            calfname_list = get_caldb_fname(caldbidxfname, telescop, instrume, detnam, filterstr, codenam)
            ncalfiles = len(calfname_list) 
            if len(calfname_list)==0 :
                self.vigfname = None
            else :
                #self.vigfname = caldbdir+os.sep+calfname_list[0]
                self.vigfname = caldbdir+os.sep+calfname_list[-1]
        elif vigfname=="NONE" :
            self.vigfname = None
        else :
            self.vigfname = vigfname

        if self.vigfname!=None :
            print("Vignet file = %s"%(self.vigfname))
            self.hzxvig = HZXVign(self.vigfname)

        return

    
if __name__=="__main__" :

    detid = 1
    detnam = "XM%02d"%(detid)

    teldefname = "../refdata/teldef/hiz_xm%02d_teldef_20251216v002.fits"%(detid)
    rmffname   = "../refdata/hzx_erosita.rmf"
    arffname   = "../refdata/hzxarf_withframe_rn0.0nm_imgall.arf"
    vigfname   = "../refdata/hzx_vig_withframe_imgall.fits"
    
    hzxcaldb =  HzxCaldb(1, teldefname, rmffname, arffname, vigfname)
    

