include ../../make.inc

testsrc = ../../../unittests
excsrc  = ../../../src
includes = -I. -I../../../src/FoX/dom/ -I../../../src/FoX/fsys/
objects = modfvsystem_test.o hmlaan_test.o olpaan_test.o hmlalon_test.o hmllolon_test.o hmlistln_test.o\
          fruit_util.o fruit.o test_helpers.o testsuite.o
cmd = $(F90)
opts = $(F90_OPTS_UT_SERIAL)

foxlib = -L../../../src/FoX/objs/lib/ -lFoX_dom -lFoX_sax -lFoX_utils -lFoX_wcml -lFoX_wxml -lFoX_fsys -lFoX_common

LIBS_T = $(foxlib) 
export SUFFIX=utser
# main program ##############################################


all: ../make.inc ../Make.common ../libraries.inc tests_serial
	rm -f libexciting.a 
tests_serial: libexciting.a $(objects) 
	$(F90) $(F90_OPTS_UT_SERIAL) -o tests_serial $(objects) libexciting_serial.a $(LIBS_T) 
        
# exciting library
libexciting.a:
	$(MAKE)  -f ../Make.common libexciting PREFIX=../../..
	ar q libexciting_serial.a `ls *.o`

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

