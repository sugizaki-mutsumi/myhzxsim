#!/usr/bin/env python

import numpy as np
from astropy.time import Time, TimeDelta


# MJDREF: 2025-01-01T00:00:00 UTC
# mjdrefi = 60676                  ### 2025-01-01
# mjdreff = 0.0008007407377590425  ### UTC -> TT time at 2025-01-01
#
# equivalent to
# reftime = Time('2025-01-01T00:00:00').tt.mjd

class MissionTime :
    def __init__(self, reftime = Time('2025-01-01T00:00:00', scale='utc')) :
        self.reftime = reftime
        self.mjdref  = reftime.tt.mjd
        self.mjdreff, self.mjdrefi = np.modf(self.mjdref)

    def mstime2asptime(self, mstimesec) :
        asptime = self.reftime + TimeDelta(mstimesec, format='sec')
        return asptime

    def asptime2mstime(self, asptime) :
        mstimesec = (asptime - self.reftime).sec
        return mstimesec
        #mstime.format='sec'


    def mjdutc2mstime(self, mjdutc) :
        asptime = Time(mjdutc, scale='utc', format='mjd')
        return self.asptime2mstime(asptime)

    def mstime2mjdutc(self, mstimesec) :
        asptime = self.mstime2asptime(mstimesec)
        return asptime.utc.mjd


    ### convert via self.mjdref
    def mstime2mjdtt(self, mstimesec) :
        mjdtt = self.mjdref + mstimesec/86400.
        return mjdtt

    def mjdtt2mstime(self, mjdtt) :
        mstimesec = (mjdtt-self.mjdref)*86400.
        return mstimesec

    ### convert via astropy Time
    def mstime2mjdtt2(self, mstimesec) :
        asptime = self.mstime2asptime(mstimesec)
        return asptime.tt.mjd

    def mjdtt2mstime2(self, mjdtt) :
        asptime = Time(mjdtt, scale='tt', format='mjd')
        return self.asptime2mstime(asptime)


### create default MissionTime as MSTime 
MSTime = MissionTime()
mstime_start = 0.0
mstime_stop  = 86400*365.25*20  ## 20 years 

if __name__=="__main__" :
    
    ### Test time
    time1 = Time('2045-01-01T00:00:00', scale='utc')

    print("MJD (UTC) = %.15lf"%(time1.utc.mjd))
    print("Time(UTC) = %s"%(time1.utc.fits))

    #print("MJD (TT)  = %.15lf"%(time1.tt.mjd))
    #print("Time(TT)  = %s"%(time1.tt.fits))

    #mjdtt = utctime.tt.mjd
    
    mstime1 = MSTime.mjdutc2mstime(time1.utc.mjd)
    print("mjdutc2mstime = %lf"%(mstime1))

    mstime2 = MSTime.mjdtt2mstime(time1.tt.mjd)
    print("mjdtt2mstime  = %lf"%(mstime2))

    mstime2 = MSTime.mjdtt2mstime2(time1.tt.mjd)
    print("mjdtt2mstime2  = %lf"%(mstime2))
    

