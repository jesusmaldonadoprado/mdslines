# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Makefile for MSDLINES
# --------------------------------------------------------------------- 
# jmaldonado at inaf-oapa // vie feb 27 14:10:56 CET 2015
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# List of objects
OBJECTS =  msdlgmr.f90 msdstypes.f90 msdmetal.f90 msdteff.f90 msdlines.f90

# Fortran compiler and path to FISIO
FC = gfortran
FITSIOLIB = /usr/local/cfitsio/lib/
FITSIOINC = /usr/local/cfitsio/include/

# Linking sequence
all: MSDLINES ;
MSDLINES:  $(OBJECTS);
	$(FC) -o MSDLINES $(OBJECTS) -L$(FITSIOLIB) -I$(FITSIOINC) -lcfitsio -lm

