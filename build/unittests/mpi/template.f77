include ../../make.inc


FC = $(F77)
FFLAGS = $(F77_OPTS_UT_MPI) 
LD = $(FC)
LDFLAGS = $(FFLAGS) $(LIBS)
AR = ar
ARFLAGS = -rc

TMPFILES = *.mod
