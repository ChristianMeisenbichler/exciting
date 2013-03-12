include ../../make.inc
testsrc = ../../../unittests

all:  tests_serial
PREFIX=../../..
include ../../Make.common

tests_serial::$(PREFIX)/build/make.inc libs $(GeneratedFiles)  $(LIBXCMOD)
	perl $(PREFIX)/build/utilities/versionstamp.pl  
	$(PREFIX)/build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p tests_serial \
	$(SOURCEDIRS)  $(LIBXCMOD)	$(testsrc)\
	&& $(MAKE) -f Makefile.mkmf tests_serial




