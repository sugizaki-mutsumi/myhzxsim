# Description of parameter files

## 1. Observation parameter file

obsparam.yml

```
### Obseravation including Crab
- obsid: 100000            # used for output filename
  objectName: XRAYSKY      # used in output fits header
  tstart: 0.0              # start time in a mission time
  tstop:  1000.0           # stop time in a mission time
  qparam: -0.58883817, 0.00332964, -0.00195868, 0.80824173   # attitude quaternion
  detlist: range(1,17).    # list of detector id number to run simulation 1,2,...,16
  ##
  mxmoncat: True   # include MAXI catalog sources (True) or not (False)
  rxscat:   True   # include ROSAT catalog sources (True) or not (False)
  xrb: DEF         # diffuse X-ray background model (Default) or None
  nxb: DEF         # non X-ray background model (Default), user-specified filename or None


### Obseravation including Galactic center
- obsid: 100001
  objectName: XRAYSKY
  tstart: 0.0
  tstop:  1000.0
  qparam: -0.82470631,  0.00189276,  0.00301748, -0.56555001
  detlist: range(1,17)
  ##
  mxmoncat: True   
  rxscat:   True   
  xrb: DEF        
  nxb: DEF        

```

## 2. X-ray source file (optional)

mytarget_grbs.yml

```
- objectName: GRB1
  ra:  83.
  dec: 13.
  xsmodel: phabs*powerlaw
  xsparam: (0.3, 1.0, 1.0)
  lcfile: grb1_lc.fits
...
```

## 3. Coniguration parameter file (optional)

The default configuration file is $HZXSIMDIR/refdata/myhzxsim_default.yml

myhzxsim_myconf.yml 
```
###. HZXSIMDIR is replaced with the value of the environment variable $HZXSIMDIR. 

#############################
##### CALDB information #####
#############################

## TELDEF
teldefdir: HZXSIMDIR/refdata/teldef
teldefver: 20251216v002

### Ancilary response (effective area) 
#arffname: HZXSIMDIR/refdata/hzxarf_withframe_rn0.0nm_imgall.arf  ### surface roughness 0 nm
arffname: HZXSIMDIR/refdata/hzxarf_withframe_rn2.0nm_imgall.arf  ### surface roughness 2 nm

### Response matrix
rmffname: HZXSIMDIR/refdata/hzx_erosita.rmf

### PSF sigma in mm
#psf_sigma: 0.2  ### = 2.29 arcmin (sigma) = 5.38 arcmin (FWHM)
psf_sigma: 0.4  ### = 4.58 arcmin (sigma) = 10.8 arcmin (FWHM) 

### Non X-ray background spectrum
nxbspecfname: HZXSIMDIR/refdata/fake_bkg_exp1e6s.pha 

### X-ray sky background spectral data in healpix image
xrbdatfname: HZXSIMDIR/simdb/xrbkg_rayspec_erobin.fits

################################
##### X-ray source catalog #####
################################
### MAXI catalog file
catmxmonfname: HZXSIMDIR/simdb/catmxmon.fits
catmxmonlcdir: HZXSIMDIR/simdb/mxgsclc

### ROSAT catalog file
catrxs2fname: HZXSIMDIR/simdb/cat2rxs_plparcut.fits
```
