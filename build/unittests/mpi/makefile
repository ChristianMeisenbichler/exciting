include ../../make.inc

testsrc = ../../../unittests
excsrc  = ../../../src
includes = -I ../../mpi/
objects = modfvsystem_test.o hmlaan_test.o \
          fruit_util.o fruit.o test_helpers.o testsuite.o
cmd = $(MPIF90)
opts = $(MPIF90_OPTS)

foxlib = -L../../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 

# main program ##############################################

all: ../make.inc ../Make.common ../libraries.inc tests_mpi
	rm -f libexciting.a
tests_mpi: libexciting.a $(objects)
	$(MPIF90) $(F90_OPTS_UT_MPI) -o tests_mpi $(objects) libexciting_mpi.a $(LIBS_T) 
        
# exciting library
libexciting.a:
	$(MAKE)  -f ../Make.common libexciting PREFIX=../../..
	ar q libexciting_mpi.a `ls *.o`
# unit tests ##############################################
	
hmlaan_test.o: $(testsrc)/hmlaan_test.F90 test_helpers.o fruit.o
	$(MPIF90) $(F90_OPTS_UT_MPI) $(includes) -c $(testsrc)/hmlaan_test.F90

modfvsystem_test.o: $(testsrc)/modfvsystem_test.F90 test_helpers.o fruit.o
	$(MPIF90) $(F90_OPTS_UT_MPI) $(includes) -c $(testsrc)/modfvsystem_test.F90
	
fruit_util.o: $(testsrc)/fruit_util.f90
	$(MPIF90) $(F90_OPTS_UT_MPI) -c $(testsrc)/fruit_util.f90

fruit.o: $(testsrc)/fruit.F90 fruit_util.o
	$(MPIF90) $(F90_OPTS_UT_MPI) -c $(testsrc)/fruit.F90

test_helpers.o: $(testsrc)/test_helpers.F90
	$(MPIF90) $(F90_OPTS_UT_MPI) $(includes) -c $(testsrc)/test_helpers.F90

testsuite.o: $(testsrc)/testsuite.F90 modfvsystem_test.o hmlaan_test.o fruit.o 
	$(MPIF90) $(F90_OPTS_UT_MPI) $(includes) -c $(testsrc)/testsuite.F90
	

# general stuff ##############################################


clean:
	rm -rf *.o *.mod tests* *.a finclude lib Makefile.* 

make.inc:
	cp ../make.inc .
Make.common:
	cp ../Make.common .
libraries.inc:
	cp ../libraries.inc .

