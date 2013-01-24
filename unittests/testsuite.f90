
program testsuite
    use fruit
    use modfvsystem_test
    use modmpi
   
    call initMPI

    CALL init_fruit

    ! matrix_tools_tests
    CALL set_unit_name ('modfvsystem - creating new matrix on 1 proc')
    CALL testNewSystem1proc
    CALL set_unit_name ('modfvsystem - creating new matrix on 4 procs')
    CALL testNewSystem4proc

    CALL fruit_summary_MPI(MPI_COMM_WORLD)

    call finitMPI

end program testsuite
