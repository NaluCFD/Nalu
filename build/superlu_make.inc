# Modified by Marc Henry de Frahan to compile SuperLU with the Intel
# MKL libraries. Also modified for build with Trilinos and Nalu.

############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _linux

#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= $(HOME)/nalu-build-intel
SUPERLULIB   	= $(SuperLUroot)/lib/libsuperlu_4.3.a

#BLASLIB   	= $(SuperLUroot)/lib/libblas.a

## ATLAS BLAS causes single-precision to fail
#BLASDEF = -DUSE_VENDOR_BLAS
#BLASLIB   	= -L/usr/lib/atlas -lblas

## This BLAS causes single-precision to fail the test in SuperLU
BLASDEF 	= -DUSE_VENDOR_BLAS
BLASLIB 	= -L/usr/lib64 -lblas

## MKL BLAS. This BLAS causes single-precision to fail the test in SuperLU
# BLASDEF 	= -DUSE_VENDOR_BLAS
# BLASLIB 	= -L${MKLROOT}/lib/intel64 -lmkl_rt
# BLASINC		= -I${MKLROOT}/include


TMGLIB       	= libtmglib.a
LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = mpicc
CFLAGS       = -O3
NOOPTS       = 
FORTRAN	     = mpif77
FFLAGS       = -O2
LOADER       = $(CC)
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#
#MATLAB	     = /usr/sww/matlab
