include ../make.inc

testsrc = ../../unittests
excsrc  = ../../src
includes = -I ../mpi/
objects = modfvsystem_test.o \
          fruit_util.o fruit.o testsuite.o 

foxlib = -L../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 

# main program ##############################################

all: tests_mpi tests_serial

tests_mpi: 
	cd mpi && make
    
tests_serial: 
	cd serial && make
    


# general stuff ##############################################

lib: libexciting_mpi.a libexciting_serial.a

libexciting_mpi.a: 
	cd ../mpi && make && ar q ../unittests/mpi/libexciting_mpi.a `ls *.o`

libexciting_serial.a: 
	cd ../serial && make && ar q ../unittests/serial/libexciting_serial.a `ls *.o`

run:
	./serial/tests_serial && mpirun -n 4 ./mpi/tests_mpi

clean:
	rm *.o *.mod tests* 

