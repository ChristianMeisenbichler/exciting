
program testsuite
    use fruit
    use modfvsystem_test
    use modHmlaan_test
    use modOlpaan_test
    use modHmlalon_test
#ifdef MPI
    use modmpi

    call initMPI
#endif
   
    CALL init_fruit

#ifdef MPI
    Call runtestcase(testcaseSystem1Proc, 'modfvsystem - System construction on 1 proc')
    Call runtestcase(testcaseSystem4Proc, 'modfvsystem - System construction on 4 procs')

    Call runtestcase(testcaseHermitianMatrixMatrix1Proc, 'modfvsystem - Hermitian matrix-matrix-multiplication on 1 proc')
    Call runtestcase(testcaseHermitianMatrixMatrix4Proc, 'modfvsystem - Hermitian matrix-matrix-multiplication on 4 procs')

    Call runtestcase(testcaseHmlaan1Proc, 'modfvsystem - Hamiltonian matrix setup (hmlaan) on 1 proc')
    Call runtestcase(testcaseHmlaan4Proc, 'modfvsystem - Hamiltonian matrix setup (hmlaan) on 4 pros')

    Call runtestcase(testcaseOlpaan1Proc, 'modfvsystem - Hamiltonian matrix setup (olpaan) on 1 proc')
    Call runtestcase(testcaseOlpaan4Proc, 'modfvsystem - Hamiltonian matrix setup (olpaan) on 4 pros')

    Call runtestcase(testcaseHmlalon1Proc, 'modfvsystem - Hamiltonian matrix setup (hmlalon) on 1 proc')
    Call runtestcase(testcaseHmlalon4Proc, 'modfvsystem - Hamiltonian matrix setup (hmlalon) on 4 pros')
#else
    Call runtestcase(testcaseSystemSerial, 'modfvsystem - System construction serial')
    Call runtestcase(testcaseHermitianMatrixMatrixSerial, 'modfvsystem - Hermitian matrix-matrix-multiplication serial')
    Call runtestcase(testcaseHmlaanSerial, 'modfvsystem - Hamiltonian matrix setup (hmlaan) serial')
    Call runtestcase(testcaseOlpaanSerial, 'modfvsystem - Overlap matrix setup (olpaan) serial')
    Call runtestcase(testcaseHmlalonSerial, 'modfvsystem - Hamiltonian matrix setup (hmlalon) serial')
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
