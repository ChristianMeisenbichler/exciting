
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
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 1 proc - is it stable?')
    CALL testHermitianMatrixMatrix1Proc_AxB_stability
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxIpC on 1 proc')
    CALL testHermitianMatrixMatrix1Proc_AxIpC
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_AxI
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_AxB
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB on 4 procs - is it stable?')
    CALL testHermitianMatrixMatrix4Proc_AxB_stability
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxIpC on 4 procs')
    CALL testHermitianMatrixMatrix4Proc_AxIpC
#else
    CALL set_unit_name ('modfvsystem - creating new complex matrix serial')
    CALL testNewComplexMatrixSerial

    CALL set_unit_name ('modfvsystem - creating new system serial')
    CALL testNewSystemSerial

    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI serial')
    CALL testHermitianMatrixMatrixSerial_AxI
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB serial')
    CALL testHermitianMatrixMatrixSerial_AxB
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxB serial - is it stable?')
    CALL testHermitianMatrixMatrixSerial_AxB_stability
    CALL set_unit_name ('modfvsystem - hermitian matrix-matrix-multiplication AxI+C serial')
    CALL testHermitianMatrixMatrixSerial_AxIpC
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
