.NOTPARALLEL:
include ../make.inc




# main program ##############################################


all:  tests_mpi tests_serial

tests_mpi:: 
	cd mpi && $(MAKE)
    
tests_serial::
	cd serial && $(MAKE)
    


# general stuff ##############################################

run: tests_mpi tests_serial
	./serial/tests_serial && mpirun -n 4 ./mpi/tests_mpi

clean:
	cd serial && $(MAKE) clean 
	cd mpi && $(MAKE) clean
#	rm serial/*.o serial/*.mod serial/*.a mpi/*.o mpi/*.mod mpi/*.a tests*

