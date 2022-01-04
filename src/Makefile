#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Note that these commands are not safe for parallel builds.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is held in BLASLIB, which is set in the
#  top-level make.inc
#
#  To remove the object files after the library is created, enter
#       make cleanobj
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov, Last updated September 30, 1994
#  ejr, May 2006.
#
#######################################################################

TOPSRCDIR = ../..
#######################################################################
#include $(TOPSRCDIR)/make.inc

SHELL = /bin/sh

#  CC is the C compiler, normally invoked with options CFLAGS.
#
CC = gcc
CFLAGS = -O3

#  Modify the FC and FFLAGS definitions to the desired compiler
#  and desired compiler options for your machine.  NOOPT refers to
#  the compiler options desired when NO OPTIMIZATION is selected.
#
#  Note: During a regular execution, LAPACK might create NaN and Inf
#  and handle these quantities appropriately. As a consequence, one
#  should not compile LAPACK with flags such as -ffpe-trap=overflow.
#
FC = gfortran
FFLAGS = -O2 -g -fbacktrace -Wall
FFLAGS_DRV = $(FFLAGS)
FFLAGS_NOOPT = -O0

#  Define LDFLAGS to the desired linker options for your machine.
#
LDFLAGS =

#  The archiver and the flag(s) to use when building an archive
#  (library).  If your system has no ranlib, set RANLIB = echo.
#
AR = ar
ARFLAGS = cr
RANLIB = ranlib

#  The location of the libraries to which you will link.  (The
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
BLASLIB      = $(TOPSRCDIR)/librefblas.a

#  DOCUMENTATION DIRECTORY
# If you generate html pages (make html), documentation will be placed in $(DOCSDIR)/explore-html
# If you generate man pages (make man), documentation will be placed in $(DOCSDIR)/man
DOCSDIR       = $(TOPSRCDIR)/DOCS
#######################################################################

.PHONY: all
all: $(BLASLIB)

ALLBLAS = M_blas.o compatible.o

$(ALLBLAS): $(FRC)

ALLOBJ = $(ALLBLAS)

$(BLASLIB): $(ALLOBJ)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

.PHONY: single double complex complex16
single: $(ALLBLAS) 
	$(AR) $(ARFLAGS) $(BLASLIB) $^
	$(RANLIB) $(BLASLIB)

double: $(ALLBLAS) 
	$(AR) $(ARFLAGS) $(BLASLIB) $^
	$(RANLIB) $(BLASLIB)

complex: $(ALLBLAS) 
	$(AR) $(ARFLAGS) $(BLASLIB) $^
	$(RANLIB) $(BLASLIB)

complex16: $(ALLBLAS) 
	$(AR) $(ARFLAGS) $(BLASLIB) $^
	$(RANLIB) $(BLASLIB)

FRC:
	@FRC=$(FRC)

.PHONY: clean cleanobj cleanlib
clean: cleanobj cleanlib
cleanobj:
	rm -f *.o
cleanlib:
	#rm -f $(BLASLIB)  # May point to a system lib, e.g. -lblas
.SUFFIXES: $(SUFFIXES) .f90 .F90 

.f90.o:
	$(FC) $(FFLAGS) -c $<
.F90.o:
	$(FC) $(FFLAGS) -c $<
