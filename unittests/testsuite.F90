
program testsuite
    use fruit
    use modfvsystem_test
    use modHmlaan_test
    use modOlpaan_test
    use modHmlalon_test
    use modOlpalon_test
    use modHmllolon_test
    use modOlplolon_test
    use modHmlistln_test
    use modOlpistln_test
    use modHamiltonandoverlapsetup_test
#ifdef MPI
    use modmpi

    call initMPI
!    Call mpi_init (MPIglobal%ierr)
#endif
   
    CALL init_fruit

#ifdef MPI
    Call runtestcase(testcaseSystem1Proc, 'modfvsystem - System construction on 1 proc')
    Call runtestcase(testcaseSystem4Proc, 'modfvsystem - System construction on 4 procs')

    Call runtestcase(testcaseHermitianMatrixMatrix1Proc, 'modfvsystem - Hermitian matrix-matrix-multiplication on 1 proc')
    Call runtestcase(testcaseHermitianMatrixMatrix4Proc, 'modfvsystem - Hermitian matrix-matrix-multiplication on 4 procs')

    Call runtestcase(testcaseHermitianmatrixIndexedUpdate1Proc, 'modfvsystem - Hermitian matrix indexed update on 1 proc')
    Call runtestcase(testcaseHermitianmatrixIndexedUpdate4Proc, 'modfvsystem - Hermitian matrix indexed update on 4 procs')

    Call runtestcase(testcaseRedistributeHermitianMatrix4Proc, 'modfvsystem - Redistribution of Hermitian Matrix on 4 procs')

    Call runtestcase(testcaseHmlaan1Proc, 'modfvsystem - Hamiltonian matrix setup (the APW-APW part / hmlaan) on 1 proc')
    Call runtestcase(testcaseHmlaan4Proc, 'modfvsystem - Hamiltonian matrix setup (the APW-APW part / hmlaan) on 4 procs')

    Call runtestcase(testcaseOlpaan1Proc, 'modfvsystem - Overlap matrix setup (the APW-APW part / olpaan) on 1 proc')
    Call runtestcase(testcaseOlpaan4Proc, 'modfvsystem - Overlap matrix setup (the APW-APW part / olpaan) on 4 procs')

    Call runtestcase(testcaseHmlalon1Proc, 'modfvsystem - Hamiltonian matrix setup (the APW-LO part / hmlalon) on 1 proc')
    Call runtestcase(testcaseHmlalon4Proc, 'modfvsystem - Hamiltonian matrix setup (the APW-LO part / hmlalon) on 4 procs')

    Call runtestcase(testcaseOlpalon1Proc, 'modfvsystem - Overlap matrix setup (the APW-LO part / olpalon) on 1 proc')
    Call runtestcase(testcaseOlpalon4Proc, 'modfvsystem - Overlap matrix setup (the APW-LO part / olpalon) on 4 procs')

    Call runtestcase(testcaseHmllolon1Proc, 'modfvsystem - Hamiltonian matrix setup (the LO-LO part / hmllolon) on 1 proc')
    Call runtestcase(testcaseHmllolon4Proc, 'modfvsystem - Hamiltonian matrix setup (the LO-LO part / hmllolon) on 4 procs')

    Call runtestcase(testcaseOlplolon1Proc, 'modfvsystem - Overlap matrix setup (the LO-LO part / olplolon) on 1 proc')
    Call runtestcase(testcaseOlplolon4Proc, 'modfvsystem - Overlap matrix setup (the LO-LO part / olplolon) on 4 procs')

    Call runtestcase(testcaseHmlistln1Proc, 'modfvsystem - Hamiltonian matrix setup (the PW-PW part / hmlistln) on 1 proc')
    Call runtestcase(testcaseHmlistln4Proc, 'modfvsystem - Hamiltonian matrix setup (the PW-PW part / hmlistln) on 4 procs')

    Call runtestcase(testcaseOlpistln1Proc, 'modfvsystem - Hamiltonian matrix setup (the PW-PW part / olpistln) on 1 proc')
    Call runtestcase(testcaseOlpistln4Proc, 'modfvsystem - Hamiltonian matrix setup (the PW-PW part / olpistln) on 4 procs')

    Call runtestcase(testcaseHamiltonandoverlapsetup1Proc, 'hamiltonandoverlapsetup on 1 proc')
    Call runtestcase(testcaseHamiltonandoverlapsetup4Proc, 'hamiltonandoverlapsetup on 4 procs')

#else
    Call runtestcase(testcaseSystemSerial, 'modfvsystem - System construction serial')
    Call runtestcase(testcaseHermitianMatrixMatrixSerial, 'modfvsystem - Hermitian matrix-matrix-multiplication serial')
    Call runtestcase(testcaseHermitianmatrixIndexedUpdateSerial, 'modfvsystem - Hermitian matrix indexed update serial')

    Call runtestcase(testcaseHmlaanSerial, 'modfvsystem - Hamiltonian matrix setup (the APW-APW part / hmlaan) serial')
    Call runtestcase(testcaseOlpaanSerial, 'modfvsystem - Overlap matrix setup (the APW-APW part / olpaan) serial')
    Call runtestcase(testcaseHmlalonSerial, 'modfvsystem - Hamiltonian matrix setup (the APW-LO part / hmlalon) serial')
    Call runtestcase(testcaseOlpalonSerial, 'modfvsystem - Overlap matrix setup (the APW-LO part / olpalon) serial')
    Call runtestcase(testcaseHmllolonSerial, 'modfvsystem - Hamiltonian matrix setup (the LO-LO part / hmllolon) serial')
    Call runtestcase(testcaseOlplolonSerial, 'modfvsystem - Overlap matrix setup (the LO-LO part / olplolon) serial')
    Call runtestcase(testcaseHmlistlnSerial, 'modfvsystem - Hamiltonian matrix setup (the PW-PW part / hmlistln) serial')
    Call runtestcase(testcaseOlpistlnSerial, 'modfvsystem - Overlap matrix setup (the PW-PW part / olpistln) serial')

    Call runtestcase(testcaseHamiltonandoverlapsetupSerial, 'hamiltonandoverlapsetup serial')

    
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    Call MPI_Finalize(MPIglobal%ierr)
!    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
