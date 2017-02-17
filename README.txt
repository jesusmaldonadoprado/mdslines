#######################################################
# README file for MSDLINES 
# version 02: vie mar 20 15:05:17 CET 2015
# jmaldonado at inaf-oapa 
#-----------------------------------------------------
# See Maldonado et al. (2015, submitted) for details.
# Please report any bug or suggestion to:
# jmaldonado@astropa.inaf.it
#######################################################

O. System requirements:
------------------------------------------------------
 It requires the gfortran compiler and the CFITSIO
 library (http://heasarc.nasa.gov/fitsio/fitsio.html).

** NOTE **********************************************
 This code was developed for HARPS / HARPS-N spectra.
 We caution that our methods are untested or unreliable
 on other spectrographs.
******************************************************* 

1. Files description:   
------------------------------------------------------ 

* 'harpn_all_4424_identified_features.dat': It contains
   the original list of identified features.

* 'stars_list.txt': File containing the list of stars
   to be calibrated. 

* 'msdlines.f90': Main program. It computes the pseudo
   equivalent widths and call the subroutines
   for stellar parameters determination.

* 'coefficients_Mstars_ratios_teff_ver02.dat':
   It contains the coefficients for computing the 
   effective temperature.  

* 'msdteff.f90': Subroutine
   to compute the effective temperature.

* 'coefficients_Mstars_ratios_sptype_ver02.dat':
   It contains the coefficients for computing the 
   spectral types.

* 'msdstypes.f90': Subroutine
   to determine the spectral type.

* 'coefficients_Mstars_ratios_metallicity_ver02.dat':
   It contains the coefficients for computing the stellar
   metallicity.

* 'msdmetal.f90': Subroutine
   to determine the stellar metallicity.

* 'msdlgmr.f90:': Subroutine for gravity, mass,
   radii, and luminosity determination. 

* 'Makefile': Makefile to compile the code

* 'gj15A.asc' and 'gj205.asc': Spectra of two stars
   to demonstrate how the code works. 

* 'README.txt': This file 
 

2. How to run the code:
------------------------------------------------------
   The list of stars to be calibrated should be included
in the 'stars_list.txt' file. The first raw of the file
should not be deleted. The second line gives the format
of the spectra. It could be either ASCII ('ascii') or
one dimensional FITS spectra ('fits'). The third raw of
the file should not be deleted. The rest of the 
lines should contain two arguments for each star:
i) a star identifier and ii) its corresponding
stellar spectra. The spectra should be in the same directory
where all the programs are. 

   The Makefile provided with the code can be used to compute
it, generating an executable file called MSDLINES (just write
make). Please, check your path to the cfitsio library and your
FORTRAN compiler. The cfitsio library is needed to read the input spectra in 
fits format. See, http://heasarc.nasa.gov/fitsio/fitsio.html
for details on how to install it if needed. 


3. Outputs:
------------------------------------------------------  
  The code produces an individual file for each star
in the 'star_list.txt' file with the extension *psdEWs.dat.
These files contain the measured pseudo equivalent widths.
It also produces three auxiliary files: 'derived_metallicities.dat'                                 
'derived_sptypes.dat', and  'derived_temperatures.dat' whose
contents are self-explanatory. Do not change the names of these
files in the codes. 

  The main output of the code is the file 'msd_summary_stellar_parameters.txt'
which contains all the derived stellar parameters. A value of "-99.0" is 
given if the star falls out of the range of applicability of our
calibrations.  

