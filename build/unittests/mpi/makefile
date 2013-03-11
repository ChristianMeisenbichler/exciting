include ../../make.inc
testsrc = ../../../unittests

all:  tests_mpi
PREFIX=../../..
include ../../Make.common

tests_mpi::$(PREFIX)/build/make.inc libs $(GeneratedFiles)  $(LIBXCMOD)
	perl $(PREFIX)/build/utilities/versionstamp.pl  
	$(PREFIX)/build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p tests_mpi \
	$(SOURCEDIRS)  $(LIBXCMOD)	$(testsrc)\
	&& $(MAKE) -f Makefile.mkmf tests_mpi




