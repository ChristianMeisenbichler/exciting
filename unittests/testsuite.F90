
program testsuite
    use fruit
    use modfvsystem_test
    use modHmlaan_test
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
#else
    Call runtestcase(testcaseSystemSerial, 'modfvsystem - System construction serial')
    Call runtestcase(testcaseHermitianMatrixMatrixSerial, 'modfvsystem - Hermitian matrix-matrix-multiplication serial')
    Call runtestcase(testcaseHmlaanSerial, 'modfvsystem - Hamiltonian matrix setup (hmlaan) serial')
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
