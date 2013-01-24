module modfvsystem_test
    use fruit
    use modfvsystem
    use modmpi

    implicit none

    real(8), parameter :: tol  = 1E-9
    real(8), parameter :: zero = 0.0D0

    contains

!------------------------------------------------------------------------------
! test testNewmatrixSerial
!------------------------------------------------------------------------------
    subroutine testNewsystem1proc
    implicit none

    integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, i, &
               world_group, tprocs_group, tprocs_comm, ierror_t, &
               proc_row, proc_col, src_rows, src_cols

    integer, parameter :: nmatp = 9
    Type (evsystem) :: system
    
    n_proc_rows_test = 1
    n_proc_cols_test = 1
    n_procs_test = n_proc_rows_test*n_proc_cols_test
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror_t)
    CALL MPI_GROUP_INCL(world_group, n_procs_test, (/(i,i=0,n_procs_test-1)/), tprocs_group, ierror_t)
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD, tprocs_group, tprocs_comm, ierror_t);

    comm = tprocs_comm
    nprocrows = n_proc_rows_test
    nproccols = n_proc_cols_test
    blocksize = 2

    CALL BLACS_GET(0, 0, context)
    CALL BLACS_GRIDINIT(context,'r', nprocrows, nproccols)

    if (rank < n_procs_test) then
      CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

      Call newsystem (system, nmatp)

      CALL assert_equals(nmatp, system%hamilton%size, 'checking newmatrix H rank')
      CALL assert_true(.not. system%hamilton%ludecomposed, 'checking newmatrix H ludecomposed')
      CALL assert_equals(nmatp, size(system%hamilton%za,1), 'checking newmatrix H size rows')
      CALL assert_equals(nmatp, size(system%hamilton%za,2), 'checking newmatrix H size cols')
      CALL assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newmatrix H=0')

      CALL assert_equals(nmatp, system%overlap%size,  'checking newmatrix S rank')
      CALL assert_true(.not. system%overlap%ludecomposed, 'checking newmatrix S ludecomposed')
      CALL assert_equals(nmatp, size(system%overlap%za,1), 'checking newmatrix S size rows')
      CALL assert_equals(nmatp, size(system%overlap%za,2), 'checking newmatrix S size cols')
      CALL assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newmatrix S=0')
    endif

    end subroutine testNewsystem1proc

!------------------------------------------------------------------------------
! test testNewmatrixParallel
!------------------------------------------------------------------------------
    subroutine testNewsystem4proc
    implicit none

    integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, i, &
               world_group, tprocs_group, tprocs_comm, ierror_t, &
               proc_row, proc_col, src_rows, src_cols

    integer, parameter :: nmatp = 9
    Type (evsystem)    :: system

    integer            :: nrows_loc, ncols_loc

    n_proc_rows_test = 2
    n_proc_cols_test = 2
    n_procs_test = n_proc_rows_test*n_proc_cols_test
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror_t)
    CALL MPI_GROUP_INCL(world_group, n_procs_test, (/(i,i=0,n_procs_test-1)/), tprocs_group, ierror_t)
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD, tprocs_group, tprocs_comm, ierror_t);

    comm = tprocs_comm
    nprocrows = n_proc_rows_test
    nproccols = n_proc_cols_test
    blocksize = 2

    CALL BLACS_GET(0, 0, context)
    CALL BLACS_GRIDINIT(context,'r', nprocrows, nproccols)

    if (rank < n_procs_test) then
      CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

write (*,*) 'rank', rank, 'row', myprocrow, 'col', myproccol

      CALL newsystem (system, nmatp)
      
      select case (rank)
        case (0)
          nrows_loc = 5
          ncols_loc = 5
        case (1)
          nrows_loc = 5
          ncols_loc = 4
        case (2)
          nrows_loc = 4
          ncols_loc = 5
        case (3)
          nrows_loc = 4
          ncols_loc = 4
      end select

      CALL assert_equals(nmatp, system%hamilton%size, 'checking newmatrix H rank')
      CALL assert_equals(nrows_loc, system%hamilton%nrows_loc, 'checking newmatrix H nrows_loc')
      CALL assert_equals(ncols_loc, system%hamilton%ncols_loc, 'checking newmatrix H ncols_loc')
      CALL assert_true(.not. system%hamilton%ludecomposed, 'checking newmatrix H ludecomposed')
      CALL assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking newmatrix H size rows')
      CALL assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking newmatrix H size cols')
      CALL assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newmatrix H=0')

      CALL assert_equals(nmatp, system%overlap%size, 'checking newmatrix S rank')
      CALL assert_equals(nrows_loc, system%overlap%nrows_loc, 'checking newmatrix S nrows_loc')
      CALL assert_equals(ncols_loc, system%overlap%ncols_loc, 'checking newmatrix S ncols_loc')
      CALL assert_true(.not. system%overlap%ludecomposed, 'checking newmatrix S ludecomposed')
      CALL assert_equals(nrows_loc, size(system%overlap%za,1), 'checking newmatrix S size rows')
      CALL assert_equals(ncols_loc, size(system%overlap%za,2), 'checking newmatrix S size cols')
      CALL assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newmatrix S=0')

    end if

    end subroutine testNewsystem4proc



end module modfvsystem_test
