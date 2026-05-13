#!/usr/bin/env python

import os
import numpy as np

from teldef import Teldef
import telcoord


hzxsimdir=os.getenv("HZXSIMDIR") 

### Crab position
ra  = 83.633212
dec = 22.01446

### Attiude Qparams for Crab region 
qparam = (-0.58883817, 0.00332964, -0.00195868, 0.80824173)
attrot = telcoord.attquat2rot(qparam)

### detector TELDEF for unit covering Crab 
detid = 7
teldefverid = 2
teldeffname = hzxsimdir+os.sep+"refdata/teldef/hiz_xm%02d_teldef_20251216v%03d.fits"%(detid, teldefverid)
this_teldef = Teldef(teldeffname) 

detxpix, detypix = telcoord.radec2detpix(ra, dec, attrot, this_teldef)

print("qapram = ", qparam)
print("Target (RA, Dec) = (%lf, %lf)"%(ra, dec))
print("Detector#%d (DETX, DETY) = (%lf, %lf)"%(detid, detxpix, detypix))


### center ra, dec
detxpix_cen = this_teldef.detxpix_cen
detypix_cen = this_teldef.detypix_cen
ra, dec = telcoord.detpix2radec(detxpix_cen, detypix_cen, attrot, this_teldef)
print("Detector#%d center (RA, Dec) = (%lf, %lf)"%(detid, ra, dec))

