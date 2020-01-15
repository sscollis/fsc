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
FFLAGS   = -r8 -O2 -col120 -c $(DEBUG)
F90FLAGS = -O2 -c $(DEBUG)
OFLAGS   = -O2 $(DEBUG) -o $(NAME)
LIB      = /usr/people/collis/lib/bslib/bslib.a
COMP     = f90

.SUFFIXES: .f90 

OBJECTS = fsc.o numrec.o

$(NAME): $(OBJECTS)
	$(COMP) $(OFLAGS) $(OBJECTS) $(LIB)

attach: attach.o numrec.o
	$(COMP) -O2 -r8 $(DEBUG) -o attach attach.o numrec.o $(LIB)

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	f77 $(FFLAGS) $*.f

clean:
	/bin/rm *.o *.mod
