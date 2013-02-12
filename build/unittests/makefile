include ../make.inc

testsrc = ../../unittests
excsrc  = ../../src
includes = -I ../mpi/
objects = modfvsystem_test.o \
          fruit_util.o fruit.o testsuite.o 

foxlib = -L../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 

# main program ##############################################

export SUFFIX=_unittests

all: make.inc Make.common libraries.inc tests_mpi tests_serial

tests_mpi: 
	cd mpi && $(MAKE)
    
tests_serial: 
	cd serial && $(MAKE)
    


# general stuff ##############################################

run:
	./serial/tests_serial && mpirun -n 4 ./mpi/tests_mpi

clean:
	cd serial && $(MAKE) clean 
	cd mpi && $(MAKE) clean
	rm -f make.inc Make.common libraries.inc
#	rm serial/*.o serial/*.mod serial/*.a mpi/*.o mpi/*.mod mpi/*.a tests*

make.inc:
	cp ../make.inc .
Make.common:
	cp ../Make.common .
libraries.inc:
	cp ../libraries.inc .
