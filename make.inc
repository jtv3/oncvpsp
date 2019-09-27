# System-dependent makefile options for ONCVPSP
# This must be carefully edited before executing "make" in src
#
# Copyright (c) 1989-2015 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
# University
 
##### Edit the following lines to correspond to your compilers ####

F77        = gfortran
F90        = gfortran
CC         = gcc
FCCPP      = cpp

FLINKER     = $(F90)

FCCPPFLAGS = -ansi -DLIBXC_VERSION=203  #This probably should not be changed

##### Edit the following optimization flags for your system ####

FFLAGS     = -O3 -ffast-math -funroll-loops
CFLAGS     = -O3		

##### Edit the following LAPACK and BLAS library paths for your system ####

LIBS = -L/usr/local/lapack/lib -llapack -lrefblas

##### Edit the following for to use libxc if available #####

OBJS_LIBXC =	exc_libxc_stub.o

# oncvpsp is compatible with libxc
# To build oncvpsp with libxc, uncomment 3 of the following lines and edit
# the paths to point to your libxc library and include directories
# make clean in src before rebuilding after changing this

##for libxc 2.1.0 and later use
LIBS += -L/home/drh/abinit/fallbacks/exports/lib -lxcf90 -lxc

##for earlier releases use
#LIBS += -L/home/drh/abinit/fallbacks/exports/lib -lxc

FFLAGS += -I/home/drh/abinit/fallbacks/exports/include
OBJS_LIBXC =	functionals.o exc_libxc.o 

