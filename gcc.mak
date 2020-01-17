#==============================================================================
#  Makefile for fsc (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 1-29-97
#
#==============================================================================
NAME     = fsc 
DEBUG    =
OPT      = -O2
FFLAGS   = -cpp -ffixed-line-length-120 -freal-4-real-8 -fdefault-real-8 \
           -fdefault-integer-8 -std=legacy $(DEFINES) $(OPT) $(DEBUG)
F90FLAGS = -cpp -freal-4-real-8 -fdefault-real-8 -fdefault-integer-8 \
           $(DEFINES) $(OPT) $(DEBUG)
OFLAGS   = $(OPT) $(DEBUG) -o $(NAME)
LIB      =
FC       = gfortran
F77      = gfortran

.SUFFIXES: .f90 

OBJECTS = fsc.o numrec.o bslib1.o bslib2.o util.o

ATTACH = attach.o numrec.o bslib1.o bslib2.o

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
	/bin/rm -f *.o *.mod attach fsc
