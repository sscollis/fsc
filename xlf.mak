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
FFLAGS   = -c -qrealsize=8 -qfixed=120 -qsuppress=cmpmsg -O2 $(DEBUG)
F90FLAGS = -c -qrealsize=8 -qsuffix=f=f90:cpp=f90 -qsuppress=cmpmsg \
-O2 -c $(DEBUG)
OFLAGS   = -O2 $(DEBUG) -o $(NAME)
LIB      = $(HOME)/lib/bslib.a
COMP     = xlf90

.SUFFIXES: .f90 

OBJECTS = fsc.o numrec.o

$(NAME): $(OBJECTS)
	$(COMP) $(OFLAGS) $(OBJECTS) $(LIB)

attach: attach.o numrec.o
	$(COMP) $(F90FLAGS) $(DEBUG) -o attach attach.o numrec.o $(LIB)

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	xlf $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod
