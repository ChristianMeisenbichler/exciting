include ../make.inc

testsrc = ../../unittests
excsrc  = ../../src
includes = -I ../mpi/
objects = modfvsystem_test.o \
          fruit_util.o fruit.o testsuite.o libexciting.a

# objects = mod_eigensystem.o modinputdom.o inputmodules.o modmpi.o timesec.o mod_timing.o modfvsystem.o \
#           modfvsystem_test.o \
#           fruit_util.o fruit.o testsuite.o 

foxlib = -L../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) $(LIB_LPK)

# main program ##############################################

tests: $(objects) 
	$(MPIF90) $(MPIF90_OPTS) -o tests $(objects) $(LIBS_T)
    


# unit tests ##############################################
	
modfvsystem_test.o: $(testsrc)/modfvsystem_test.f90 fruit.o
	$(MPIF90) $(MPIF90_OPTS) $(includes) -c $(testsrc)/modfvsystem_test.f90
	
fruit_util.o: $(testsrc)/fruit_util.f90
	$(MPIF90) $(MPIF90_OPTS) -c $(testsrc)/fruit_util.f90

fruit.o: $(testsrc)/fruit.f90 fruit_util.o
	$(MPIF90) $(MPIF90_OPTS) -c $(testsrc)/fruit.f90

testsuite.o: $(testsrc)/testsuite.f90 modfvsystem_test.o fruit.o 
	$(MPIF90) $(MPIF90_OPTS) $(includes) -c $(testsrc)/testsuite.f90
	

# general stuff ##############################################

lib: libexciting.a

libexciting.a: 
	cd ../mpi && make && ar q ../unittests/libexciting.a `ls *.o`

run:
	mpirun -n 4 ./tests

clean:
	rm *.o *.mod tests libexciting.a

