#==============================================================================
#  Makefile for fsc (Cray)
#
#  Author:  Scott Collis
#
#  Revised: 1-29-97
#
#==============================================================================
NPROC    = 4
NAME     = fsc 
DEBUG    = 
FFLAGS   = -N80 -c $(DEBUG)
F90FLAGS = -c $(DEBUG)
OFLAGS   = $(DEBUG) -o $(NAME)
LIB      = -limsl
COMP     = f90

.SUFFIXES: .f90 

OBJECTS = fsc.o numrec.o

$(NAME): $(OBJECTS)
	$(COMP) $(OFLAGS) $(OBJECTS) $(LIB)

.f90.o:
	$(COMP) $(F90FLAGS) $*.f90 

.f.o:
	f77 $(FFLAGS) $*.f

clean:
	/bin/rm *.o
