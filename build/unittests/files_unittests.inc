hmlistln_test.o: $(testsrc)/hmlistln_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/hmlistln_test.F90

hmllolon_test.o: $(testsrc)/hmllolon_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/hmllolon_test.F90

hmlalon_test.o: $(testsrc)/hmlalon_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/hmlalon_test.F90

olpaan_test.o: $(testsrc)/olpaan_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/olpaan_test.F90

hmlaan_test.o: $(testsrc)/hmlaan_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/hmlaan_test.F90

modfvsystem_test.o: $(testsrc)/modfvsystem_test.F90 test_helpers.o fruit.o
	$(cmd) $(opts) $(includes) -c $(testsrc)/modfvsystem_test.F90
	
fruit_util.o: $(testsrc)/fruit_util.f90
	$(cmd) $(opts) -c $(testsrc)/fruit_util.f90

fruit.o: $(testsrc)/fruit.F90 fruit_util.o
	$(cmd) $(opts) -c $(testsrc)/fruit.F90

test_helpers.o: $(testsrc)/test_helpers.F90
	$(cmd) $(opts) $(includes) -c $(testsrc)/test_helpers.F90

testsuite.o: $(testsrc)/testsuite.F90 modfvsystem_test.o hmlaan_test.o fruit.o 
	$(cmd) $(opts) $(includes) -c $(testsrc)/testsuite.F90
