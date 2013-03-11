export SUFFIX=mpi
include ../make.inc

all:
ifeq ($(BUILDMPI),true)
	$(MAKE)  -f ../Make.common PREFIX=../..
endif

libs:
ifeq ($(BUILDMPI),true)
	$(MAKE) libs -f ../Make.common PREFIX=../..

endif

clean:
	$(MAKE) clean -f ../Make.common PREFIX=../..

cleanlibs:
	$(MAKE) cleanlibs  -f ../Make.common PREFIX=../..





	
