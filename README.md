# myhzxsim

<!--
## *Cautions: Descriptions below are mostly under construction*
<br>
-->

myhzxsim is a simulator of X-ray astronomical observation with wide-field X-ray telescope utilizing lobster-eye micropore optics (MPO). Its primary purpose is research and development for X-ray satellite missions, particularly the HiZ-GUNDAM mission.
The software is written in Python and works with the HEASOFT standard high-energy astronomical library. All of the codes are tested on Mac and Linux (Ubuntu) using the Anaconda Python environment.



## Requirements

- Python(tested in version 3.12)
  - numpy, scipy, astropy, matplotlib, ...
  - healpy (https://healpy.readthedocs.io/en/latest/)
  - regions (https://astropy-regions.readthedocs.io/en/stable/)
- HEASOFT including PyXspec (tested in version 6.35, 6.36)
- (CALDB)

## How to run
### 1. [Quick start](./doc/Quickstart.md)
### 2. [Extended tutorials](./doc/Tutorials.md)
### 3. [Descriptions of input parameter files](./doc/Parameters.md)

## Major points
### 1. Calibration database (CALDB) files

  CALDB files of telescope definition (TELDEF), response matrix (RMF), ancillary response (ARF), ... are stored in the directory `refdata/`.
  In the future, these files will be moved to the location pointed by the CALDB environment variable. 

### 2. Response functions of Lobster-eye MPO

#### Point Spread Function (PSF)
  myhzxsim employs point response data of Lobster-eye MPO calculated by Geant4 X-ray tracing simulations.
  To take account of X-ray photons coming from X-ray sources located outside of the FOV, the response data is extended to the image region twice larger than the nominal FOV.
  The data is stored in a directory, `simdb/psfdb/`.<br>
  The shadows of support structures such as frames and baffles, are calculated for each source according to the incident angle.<br>
  The blurring of the PSF in real instruments due to imperfect optical alignment is implemented by a Gaussian function.

#### Effective area (ancillary response file)
  The MPO  effective area is calculated with Geant4 X-ray tracing simulation and also analytic methods for ideal devices.
  The factor of imperfectness in real instruments is parameterized by the reflectivity degradation due to the surface roughness (0 or 2 nm). 

### 3. Response functions of (CCD) X-ray imagers
The present tentative RMF is based on those for the eROSITRA pnCCD Data-Release 1 (https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_arfrmf/). 

### 4. Soft X-ray background emission
  <!-- One default option.  Select the default or none.  -->
  
  Soft X-ray backgrounds consist of unresolved extragalactic sources (mostly AGNs) and Galactic diffuse emissions that are thought to come from local hot bubbles and Galactic halos. 
  Using the [ROSAT all-sky survey data and its background tool (sxrbg) released via HEASARC](https://heasarc.gsfc.nasa.gov/Tools/xraybg_help.html), these X-ray spectra are modeled for individual 1x1 deg<sup>2</sup> sky regions pixelized by healpix of Nside = 64 (total 49152 pixels). The extragalactic component is known to be uniform and approximated by a power-law function with a photon index ~ 1.5. 
  The residual Galactic component was fitted with a thin-thermal plasma emission model.
  The obtained model data is stored in `simdb/xrbkg_rayspec.fits`. It is used in default.

### 5. Non X-ray backgound (NXB)
 
  Non X-ray background (NXB) represents residual cosmic-ray events after screening these events with event grades (pixel hit patters), pulse-height upper/lower limits.
  The default file is refdata/fake_bkg_exp1e6s.pha, which assumes a flat spectrum of ~0.005 counts cm<sup>-2</sup>  s<sup>-1</sup> keV<sup>-1</sup> throughout the energy band (e.g. Zhang+2025).

### 6. X-ray catalog sources
There are two catalog data sets. Choose both, one of them, or neither.

#### MAXI monitoring source<br>
The MAXI team releases daily light curves of ~400 bright variable / flaring sources via the [MAXI home page](http://maxi.riken.jp/). <!--These targets are selected based on daily image. -->
Among them, ~100 objects are constantly detected in a daily time scale.
They are mostly Galactic X-ray binaries embedding black-hole or neutron stars. A half of them are located on the Galactic plane.
These soft X-ray spectra suffer heavy interstellar absorption. <br>
As for these 100 objects, MAXI light curves in 2-4 keV band and spectral models obtained from past observation results are utilized. <br>
The catalog data is stored in a file, simdb/catmxmon.fits.


#### ROSAT all-sky source catalog (2RXS)<br>
2RXS catalog is the latest (last?) release of the ROSAT all-sky source catalogs, including 135,000 sources with ~30% possible spurious detections (Boller+2016). Too bright objects (Sco X-1, …) and regions (Cygnus Vela, ..) are excluded.
The data contains results of spectral fits with power-law, blackbody, and thin-thermal models to bright objects with enough photon statistics (roughly > 100 photons).<br>
In myhzxsim, power-law spectral models are utilized if these model parameters are available.
These catalog data are compiled into a file, `simdb/cat2rxs_plparcut.fits`.<br>
Some of 2RXS objects are included in the MAXI catalog. In these objects, data in the MAXI catalog is used.


### 7. Auxillary data (satellite attitude and orbit)
  To run simulation, attitude quaternion parameters are required.

  ***To be described later***

### 8. Output file
  myhzxsim generates a processed event FITS file including columns of
  TIME, PI, X, Y, RAWX, RAWY, DETX, DETY, RA, DEC, PHA, PI, RA_OBJ, DEC_OBJ, PHOTON_E.

  ***The detail will be described later***

<!-- 
| | |
|---|---|
|RA_OBJ  | Source sky position RA in degre (-999. for NXB)|
|DEC_OBJ | Source sky position Dec in degree (-999. for NXB)|
|PHOTON_E | Photon energy in keV   (-1. for NXB)|
| | |
-->