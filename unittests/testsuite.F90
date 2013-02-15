
program testsuite
    use fruit
    use modfvsystem_test
#ifdef MPI
    use modmpi
#endif
   
    CALL init_fruit

#ifdef MPI
    call initMPI
#endif

#ifdef MPI
    CALL set_unit_name ('modfvsystem - creating new complex matrix on 1 proc')
    CALL testNewComplexMatrix1proc
    CALL set_unit_name ('modfvsystem - creating new complex matrix on 4 procs')
    CALL testNewComplexMatrix4proc

    CALL set_unit_name ('modfvsystem - creating new system on 1 proc')
    CALL testNewSystem1proc
    CALL set_unit_name ('modfvsystem - creating new system on 4 procs')
    CALL testNewSystem4proc

    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI on 1 proc')
    CALL testHermitianMatrixMatrix1Proc_AxI
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 1 proc')
    CALL testHermitianMatrixMatrix1Proc_AxB
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 1 proc, square matrices')
    CALL testHermitianMatrixMatrix1Proc_AxB_square
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication IxI+C on 1 proc')
    CALL testHermitianMatrixMatrix1Proc_IxIpC
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_AxI
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_AxB
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 4 procs, square matrices')
    CALL testHermitianMatrixMatrix4Proc_AxB_square
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication IxI+C on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_IxIpC
#else
    CALL set_unit_name ('modfvsystem - creating new complex matrix serial')
    CALL testNewComplexMatrixSerial

    CALL set_unit_name ('modfvsystem - creating new system serial')
    CALL testNewSystemSerial

    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI serial')
    CALL testHermitianMatrixMatrixSerial_AxI
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB serial')
    CALL testHermitianMatrixMatrixSerial_AxB
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB serial, square general matrices')
    CALL testHermitianMatrixMatrixSerial_AxB_square
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication IxI+C serial')
    CALL testHermitianMatrixMatrixSerial_IxIpC
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
