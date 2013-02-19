include ../../make.inc

testsrc = ../../../unittests
excsrc  = ../../../src
includes = -I ../../mpi/
objects = modfvsystem_test.o hmlaan_test.o \
          fruit_util.o fruit.o test_helpers.o testsuite.o
cmd = $(MPIF90)
opts = $(F90_OPTS_UT_MPI)

foxlib = -L../../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 

# main program ##############################################

all: ../make.inc ../Make.common ../libraries.inc tests_mpi
	rm -f libexciting.a
tests_mpi: libexciting.a $(objects)
	$(cmd) $(opts) -o tests_mpi $(objects) libexciting_mpi.a $(LIBS_T) 
        
# exciting library
libexciting.a:
	$(MAKE)  -f ../Make.common libexciting PREFIX=../../..
	ar q libexciting_mpi.a `ls *.o`

# unit tests ##############################################
include ../files_unittests.inc

# general stuff ##############################################


clean:
	rm -rf *.o *.mod tests* *.a finclude lib Makefile.* 

make.inc:
	cp ../make.inc .
Make.common:
	cp ../Make.common .
libraries.inc:
	cp ../libraries.inc .

