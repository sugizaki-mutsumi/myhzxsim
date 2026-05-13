#!/usr/bin/env python

import os
import numpy as np

import astropy.io.fits as pyfits


### GRB light curve model by exponential function ###
###  norm:   count rate (c/s) at t=0
###  tau :   exponential decay time
###  tstart: start time (in obs time 0--1000 s) 
grb_norm_tau_tstart_list = [
    [10,  10,  100],
    [2,   100, 200],
    [100, 10,  300],
]


vtime = np.linspace(0, 1000, 1001)  ### time from 0 to 1000 s


num = len(grb_norm_tau_tstart_list)
for idx in range(num) :
    outfname = "grb%d_lc.fits"%(idx+1)
    if os.path.isfile(outfname) :
        print("output file %s already exists"%(outfname))
        continue
    norm   = grb_norm_tau_tstart_list[idx][0]
    tau    = grb_norm_tau_tstart_list[idx][1]
    tstart = grb_norm_tau_tstart_list[idx][2]
    
    vrate = np.exp(-1.0/tau * (vtime-tstart))*norm
    vrate *= (vtime>=tstart)

    col1 = pyfits.Column(name='TIME', format='D', unit='s',  array=vtime)
    col2 = pyfits.Column(name='RATE', format='D', unit='c/s', array=vrate)
    
    cols = pyfits.ColDefs([col1, col2])
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(outfname)
    print("output file = %s has been created"%(outfname))
    
