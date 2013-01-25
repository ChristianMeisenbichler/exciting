
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
    CALL set_unit_name ('modfvsystem - creating new matrix on 1 proc')
    CALL testNewSystem1proc
    CALL set_unit_name ('modfvsystem - creating new matrix on 4 procs')
    CALL testNewSystem4proc
#else
    CALL set_unit_name ('modfvsystem - creating new matrix serial')
    CALL testNewSystemserial
#endif

#ifdef MPI
    CALL fruit_summary_MPI(MPI_COMM_WORLD)
    call finitMPI
#else
    CALL fruit_summary
#endif


end program testsuite
