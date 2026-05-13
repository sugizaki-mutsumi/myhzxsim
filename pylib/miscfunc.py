#!/usr/bin/env python

import os, sys
import math
import numpy as np
import astropy.io.fits as pyfits
import astropy.wcs as wcs


def write_wcsfile(w, outfname, write_pole=False) :

    if os.path.isfile(outfname) :
        os.remove(outfname)
    outfile = open(outfname, "w")

    if write_pole :
        outfile.write("LATPOLE = %lf\n"%(w.wcs.latpole))
        outfile.write("LONPOLE = %lf\n"%(w.wcs.lonpole))

    outfile.write("CTYPE1 = '%s'\n" %(w.wcs.ctype[0]))
    outfile.write("CUNIT1 = '%s'\n" %(w.wcs.cunit[0]))
    outfile.write("CRPIX1 = %lf\n"%(w.wcs.crpix[0]))
    outfile.write("CRVAL1 = %lf\n"%(w.wcs.crval[0]))
    outfile.write("CDELT1 = %lf\n"%(w.wcs.cdelt[0]))

    outfile.write("CTYPE2 = '%s'\n" %(w.wcs.ctype[1]))
    outfile.write("CUNIT2 = '%s'\n" %(w.wcs.cunit[1]))
    outfile.write("CRPIX2 = %lf\n"%(w.wcs.crpix[1]))
    outfile.write("CRVAL2 = %lf\n"%(w.wcs.crval[1]))
    outfile.write("CDELT2 = %lf\n"%(w.wcs.cdelt[1]))

    outfile.write("PC1_1  = %lf\n"%(w.wcs.pc[0][0]))
    outfile.write("PC1_2  = %lf\n"%(w.wcs.pc[0][1]))
    outfile.write("PC2_1  = %lf\n"%(w.wcs.pc[1][0]))
    outfile.write("PC2_2  = %lf\n"%(w.wcs.pc[1][1]))
    outfile.close()

