#!/usr/bin/env python

import os, sys
import numpy as np

import yaml

hzxsimdir = os.getenv("HZXSIMDIR")
sys.path.append(hzxsimdir+os.sep+"pylib")
simdbdir   = hzxsimdir+os.sep+"simdb"
refdatadir = hzxsimdir+os.sep+"refdata"


from hzxsim import HZXSim

from hzxpsf import HZXPSF, get_psffilelist
from hzxcaldb import HzxCaldb
from evtproc import EventProc


from xraycatalog import CatmxmonSources, Cat2rxsSources


if __name__=="__main__" :


    if len(sys.argv)<2 :
        print("(Usage) myhzxsim.py (Obs parameter filename) [-s (X-ray source filname)] [-c (config filename)]")
        sys.exit()


        
    ### obs. parameter file
    obsparfname = sys.argv[1]
    ### read observation parameter
    with open(obsparfname) as f :
        obsparams = yaml.safe_load(f)


    ### default options
    xraysrcfname = None
    defparfname  = None

    narg = len(sys.argv)
    for idx in range(2, narg) :
        if sys.argv[idx]=="-s" :
            xraysrcfname = sys.argv[idx+1]
            idx+=1
        if sys.argv[idx]=="-c" :
            defparfname = sys.argv[idx+1]
            idx+=1


            
    ### optional source list
    xraysrclist = None ### default
    if xraysrcfname!=None :
        if os.path.isfile(xraysrcfname) :
            print("Reading xray source list %s"%(xraysrcfname))
            with open(xraysrcfname) as f :
                xraysrclist = yaml.safe_load(f)
        

    ### default configuration parameter file
    if defparfname == None :
        defparfname = hzxsimdir+os.sep+"refdata/myhzxsim_default.yml"


    with open(defparfname) as f :
        defparams = yaml.safe_load(f)
    for key in defparams.keys() :
        if isinstance(defparams[key], str) :
            if defparams[key].find("HZXSIMDIR")>=0 :
                defparams[key] = defparams[key].replace("HZXSIMDIR", hzxsimdir)
            print(key, defparams[key])

            
    
    ### MAXI catalog
    catfname = defparams["catmxmonfname"]
    catmxmonSources = CatmxmonSources(catfname)
    mxlcdir_local = defparams["catmxmonlcdir"]

    
    ### ROSAT 2RXS catalog
    catfname = defparams["catrxs2fname"]
    cat2rxsSources = Cat2rxsSources(catfname)
    
    #catfname = simdbdir+os.sep+"cat2rxs_plparcut.fits"


    ### X-ray background emission data    
    xrbdatfname_default  = defparams["xrbdatfname"]
    ##xrbdatfname_default  = simdbdir+os.sep+"xrbkg_rayspec_erobin.fits"
        
    ### Non-X-ray background spectrum    
    nxbspecfname_default = defparams["nxbspecfname"]
    #nxbspecfname_default = refdatadir+os.sep+"fake_bkg_exp1e6s.pha"

    
    
    ### Response matrix
    rmffname = defparams["rmffname"]
    #rmffname = refdatadir+os.sep+"hzx_erosita.rmf"

    ### Ancilary response (effective area) 
    arffname = defparams["arffname"]
    #arffname = refdatadir+os.sep+"hzxarf_withframe_rn0.0nm_imgall.arf"
    
    ### Vignetting data default
    vigfname   = refdatadir+os.sep+"hzx_vig_withframe_imgall.fits"


    ### Teldef
    teldefdir = defparams["teldefdir"]
    teldefver = defparams["teldefver"]
    #teldefdir = refdatadir+os.sep+"teldef"
    #teldefver = "20251216v002"

    
    ### HZX simulator and event processor
    hzxsim = HZXSim()
    evtproc = EventProc()   


    ### init PSF database
    psffitsdir =  refdatadir+os.sep+"psfdb"
    ekevlist = np.append( np.linspace(0.2, 3.0, 29), np.linspace(3.2, 8.0, 25) )
    this_elist, this_flist = get_psffilelist(ekevlist, psffitsdir)

    
    #psf_sigma = 0.2 # mm
    psf_sigma = defparams["psf_sigma"]
    hzxpsfdb = HZXPSF(sigma=psf_sigma, elist=this_elist, flist=this_flist, psfdir=psffitsdir)
    hzxsim.setHzxpsf(hzxpsfdb)


    #### quaternion parameters from attitude file
    #inv=True  ### it is turned to be befault.

    for obsparam in obsparams :
    
        obsid  = obsparam["obsid"]
        if not isinstance(obsid, str) : obsid = "%d"%(obsid)

        
        objectName = obsparam["objectName"]
        tstart = float(obsparam["tstart"])
        tstop = float(obsparam["tstop"])
    
        detid_list = obsparam["detlist"]
        if isinstance(detid_list, str) : detid_list = eval(detid_list)
        num_detid = len(detid_list)
    
        qparam = obsparam["qparam"]
        if isinstance(qparam, str) : qparam = eval(qparam)


        mxmoncat = obsparam["mxmoncat"]
        rxscat   = obsparam["rxscat"]

        xrb = obsparam["xrb"]
        nxb = obsparam["nxb"]

        
        ### Soft X-ray background
        if xrb.upper()[:3]=="DEF" :
            xrbdatfname = xrbdatfname_default
            is_xrb = True
        else :
            is_xrb = False

        #### Non X-ray background
        if nxb.upper()[:3]=="DEF" :
            nxbspecfname = nxbspecfname_default
            is_nxb = True
        elif nxb.upper()[:3]!="NON" :
            nxbspecfname = nxb
            if not os.path.isfile(nxbspecfname) :
                print("Warning: nxbspecfname = %s does not exist."%(nxbspecfname))
                is_nxb = False
        else :
            is_nxb = False
            

        if "outprefix" in obsparam.keys() :
            outprefix = obsparam["outprefix"]
        else :
            outprefix = "hzx"
            
        if "outsuffix" in obsparam.keys() :
            outsuffix = obsparam["outsuffix"]
        else :
            outsuffix = ""
            
        ### output directory name
        if "outdir" in obsparam.keys() :
            outdir = obsparam["outdir"]
        else :
            outdir = "DEF"
        if outdir == "DEF" :
            outdir = obsid

        ### 
        hzxsim.setAttQuat(qparam)

        for detid in detid_list :

            if not os.path.isdir(outdir) :
                os.mkdir(outdir)

            outfname = outdir+os.sep+"%s%sxm%02d%s.evt"%(outprefix, obsid, detid, outsuffix)
            if os.path.isfile(outfname) :
                print("Output file = %s already exists"%(outfname))
                continue
            
            print("XM02%d"%(detid))
    
            ### init CALDB
            #hzxcaldb = Hzxcaldb(cmosid)
            #hzxsim.setCaldb(hzxcaldb)
            #teldef = hzxcaldb.teldef
            #backscale = teldef.det_xsiz*teldef.det_xscl * teldef.det_ysiz*teldef.det_yscl * 0.01 # (mm^2->cm^2) 

            ### init CALDB
            teldefname = teldefdir+os.sep+"hiz_xm%02d_teldef_%s.fits"%(detid, teldefver)
            hzxcaldb =  HzxCaldb(detid, teldefname, rmffname, arffname, vigfname)
            hzxsim.setCaldb(hzxcaldb)

            
            ### init Event data
            hzxsim.initEventData()
    

            ### User X-ray source list
            if xraysrclist != None :
                hzxsim.simEventGen_XraysrcList(xraysrclist, tstart, tstop)
            
            
            ### MAXI monitor sources
            if mxmoncat :
                hzxsim.simEventGen_MaximoniSources(catmxmonSources, tstart, tstop, mxlcdir=mxlcdir_local)
    
            ### ROSAT sources
            if rxscat :
                hzxsim.simEventGen_RosatSources(cat2rxsSources, tstart, tstop, mxmonSources=catmxmonSources)
    
    
            ## X-ray background events (from ROSAT all-sky survey data) 
            if is_xrb :
                hzxsim.initXRBKSpec(xrbdatfname)
                hzxsim.simEventGen_XRBK(tstart, tstop) 
    
    
            ## Non-X-ray background events (from simulation)
            if is_nxb :
                hzxsim.readNXBKSpec(nxbspecfname, backscale=1)
                hzxsim.simEventGen_NXBK(tstart, tstop) 
    
            
            ### Detector response simulation
            hzxsim.simDetectorResp()
    
            ### Init data processing
            #hzxevtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, inv)
            evtproc.setObsParams(hzxcaldb.teldef, tstart, tstop, qparam, objectName, obsid)
            ### Store date to file
            evtproc.procEventToFile(hzxsim.vra, hzxsim.vdec, hzxsim.vtime, hzxsim.vdetxpix, hzxsim.vdetypix, hzxsim.vpha, hzxsim.vene, outfname)

