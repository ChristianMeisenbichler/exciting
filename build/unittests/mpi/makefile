include ../../make.inc

testsrc = ../../../unittests
excsrc  = ../../../src
includes = -I ../../mpi/
objects = modfvsystem_test.o \
          fruit_util.o fruit.o testsuite.o 
cmd = $(MPIF90)
opts = $(MPIF90_OPTS)

foxlib = -L../../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 

# main program ##############################################

all: tests_mpi
tests_mpi: $(objects) 
	$(cmd) $(opts) -o tests_mpi $(objects) libexciting_mpi.a $(LIBS_T) 
        

# unit tests ##############################################
	
modfvsystem_test.o: $(testsrc)/modfvsystem_test.F90 fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/modfvsystem_test.F90
	
fruit_util.o: $(testsrc)/fruit_util.f90
	$(cmd) $(opts) -c $(testsrc)/fruit_util.f90

fruit.o: $(testsrc)/fruit.F90 fruit_util.o
	$(cmd) $(opts) -c $(testsrc)/fruit.F90

testsuite.o: $(testsrc)/testsuite.F90 modfvsystem_test.o fruit.o 
	$(cmd) $(opts) $(includes) -c $(testsrc)/testsuite.F90
	

# general stuff ##############################################


clean:
	rm *.o *.mod tests* 

