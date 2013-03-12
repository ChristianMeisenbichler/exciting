include ../make.inc
export SUFFIX=smp
ifeq ($(BUILDSMP),true)
MTCOMPILE=true
endif

all:
ifeq ($(BUILDSMP),true)
	$(MAKE)  -f ../Make.common PREFIX=../..

endif
	

libs:
ifeq ($(BUILDSMP),true)
	$(MAKE) libs -f ../Make.common PREFIX=../..
endif

clean:
	$(MAKE) clean -f ../Make.common PREFIX=../..

cleanlibs:
	$(MAKE) cleanlibs  -f ../Make.common PREFIX=../..




	
