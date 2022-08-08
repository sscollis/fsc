#==============================================================================
#  Makefile for fsc (GCC) 
#
#  Author:  Scott Collis
#
#  Revised: 2-20-2020 
#
#==============================================================================
#
# If you have Numerical Recipes with a valid licence you can activate
# this. Of course, this commercial code is not distributed with FSC
# To use you must supply the needed routines in numrec.f and you must
# have a valid license to use that code.
#
# Note:  This is completely optional as inline integration routines work
#        well.  This is just provided for historical reference.
#
SHELL := /bin/bash

ifdef USE_NR
  DEFINES += -DUSE_NR
  #NROBJ = nr_runge.o
  NRLIB = -L../NR-utilities -lnr
endif
#
NAME     = fsc
DEBUG    =
OPT      = -O2
FFLAGS   = -cpp -ffixed-line-length-120 -freal-4-real-8 -fdefault-real-8 \
           -fdefault-integer-8 -std=legacy $(DEFINES) $(OPT) $(DEBUG)
F90FLAGS = -cpp -freal-4-real-8 -fdefault-real-8 -fdefault-integer-8 \
           $(DEFINES) $(OPT) $(DEBUG)
OFLAGS   = $(OPT) $(DEBUG) -o $(NAME)
LIB      = $(NRLIB)
FC       = gfortran
F77      = gfortran

.SUFFIXES: .f90 

OBJECTS = fsc.o $(NROBJ) bslib1.o bslib2.o util.o

ATTACH = attach.o $(NROBJ) bslib1.o bslib2.o util.o

all: fsc attach

$(NAME): $(OBJECTS)
	$(FC) $(OFLAGS) $(OBJECTS) $(LIB)

attach: $(ATTACH) 
	$(FC) $(OPT) $(DEBUG) $(ATTACH) $(LIB) -o attach

.f90.o:
	$(FC) $(F90FLAGS) -c $*.f90 

.f.o:
	$(F77) $(FFLAGS) -c $*.f

clean:
	$(RM) *.o *.mod fsc attach

distclean:
	make clean
	./cleanup.sh

check:
	./fsc < test.inp && \
	which ndiff && \
	./ndiff -h && \
	./ndiff --absolute-error 1e-8 --relative-error 1e-8 cprofile.dat cprofile.ref && \
	./ndiff --absolute-error 1e-8 --relative-error 1e-8 sprofile.dat sprofile.ref; \
	STATUS=$$?;\
	echo "ndiff existed with $$STATUS"; \
	exit $$STATUS
