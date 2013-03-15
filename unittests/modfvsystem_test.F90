! !MODULE:  modfvsystem_test
! !DESCRIPTION:
! Modules with test routines for basic fv system functionality
!
! !REVISION HISTORY:
!   Created: Februrary 2013 (G. Huhs - BSC)
!
module modfvsystem_test
    use fruit
    use test_helpers
    use modfvsystem
#ifdef MPI
    use modmpi
#endif

    implicit none

    contains

    subroutine testcaseSystemSerial
      Implicit None

      CALL set_test_name ('new complex matrix')
      CALL testNewComplexMatrixSerial

      CALL set_test_name ('new hermitian matrix')
      CALL testNewHermitianMatrixSerial

      CALL set_test_name ('new system')
      CALL testNewSystemSerial

    end subroutine testcaseSystemSerial


#ifdef MPI
    subroutine testcaseSystem1Proc
      Implicit None

      CALL set_test_name ('new complex matrix')
      CALL testNewComplexMatrix1proc

      CALL set_test_name ('new hermitian matrix')
      CALL testNewHermitianMatrix1proc

      CALL set_test_name ('new system')
      CALL testNewSystem1proc

    end subroutine testcaseSystem1Proc
#endif


#ifdef MPI
    subroutine testcaseSystem4Proc
      Implicit None

      CALL set_test_name ('new complex matrix')
      CALL testNewComplexMatrix4proc
!       CALL set_test_name ('new complex matrix, distribute rows')
!       CALL testNewComplexMatrix4procDistRows
      CALL set_test_name ('new complex matrix, distribute cols')
      CALL testNewComplexMatrix4procDistCols

      CALL set_test_name ('new hermitian matrix')
      CALL testNewHermitianMatrix4proc
!       CALL set_test_name ('new hermitian matrix, distribute rows')
!       CALL testNewHermitianMatrix4procDistRows
      CALL set_test_name ('new hermitian matrix, distribute cols')
      CALL testNewHermitianMatrix4procDistCols
      CALL set_test_name ('new system')
      CALL testNewSystem4proc
      CALL set_test_name ('new system, distribute cols')
      CALL testNewSystem4procDistCols

    end subroutine testcaseSystem4Proc
#endif


    subroutine testcaseHermitianMatrixMatrixSerial
      Implicit None

      CALL set_test_name ('AxI')
      CALL testHermitianMatrixMatrixSerial_AxI
      CALL set_test_name ('AxB')
      CALL testHermitianMatrixMatrixSerial_AxB
      CALL set_test_name ('AxB, square general matrices')
      CALL testHermitianMatrixMatrixSerial_AxB_square
      CALL set_test_name ('IxI+C')
      CALL testHermitianMatrixMatrixSerial_IxIpC

      CALL set_test_name ('AxB+C big hamiltonian')
      CALL testHermitianMatrixMatrixSerial_AxBpC_bighamiltonian

    end subroutine testcaseHermitianMatrixMatrixSerial


#ifdef MPI
    subroutine testcaseHermitianMatrixMatrix1Proc
      Implicit None

      CALL set_test_name ('AxI')
      CALL testHermitianMatrixMatrix1Proc_AxI
      CALL set_test_name ('AxB')
      CALL testHermitianMatrixMatrix1Proc_AxB
      CALL set_test_name ('AxB square matrices')
      CALL testHermitianMatrixMatrix1Proc_AxB_square
      CALL set_test_name ('IxI+C')
      CALL testHermitianMatrixMatrix1Proc_IxIpC

      CALL set_test_name ('AxB+C big hamiltonian')
      CALL testHermitianMatrixMatrix1Proc_AxBpC_bighamiltonian

    end subroutine testcaseHermitianMatrixMatrix1Proc
#endif


#ifdef MPI
    subroutine testcaseHermitianMatrixMatrix4Proc
      Implicit None

      CALL set_test_name ('AxI')
      CALL testHermitianMatrixMatrix4Proc_AxI
      CALL set_test_name ('AxB')
      CALL testHermitianMatrixMatrix4Proc_AxB
      CALL set_test_name ('AxB, square matrices')
      CALL testHermitianMatrixMatrix4Proc_AxB_square
      CALL set_test_name ('IxI+C')
      CALL testHermitianMatrixMatrix4Proc_IxIpC
      CALL set_test_name ('AxB+C with col distribution')
      CALL testHermitianMatrixMatrix4Proc_AxBpC

      CALL set_test_name ('AxB+C with big hamiltonian')
      CALL testHermitianMatrixMatrix4Proc_AxBpC_bighamiltonian
    end subroutine testcaseHermitianMatrixMatrix4Proc
#endif


    subroutine testcaseHermitianmatrixIndexedUpdateSerial
      Implicit None

      CALL set_test_name ('Trying to update whole matrix')
      CALL testHermitianmatrixIndexedUpdateSerial

    end subroutine testcaseHermitianmatrixIndexedUpdateSerial


#ifdef MPI
    subroutine testcaseHermitianmatrixIndexedUpdate1Proc
      Implicit None

      CALL set_test_name ('Trying to update whole matrix')
      CALL testHermitianmatrixIndexedUpdate1Proc

    end subroutine testcaseHermitianmatrixIndexedUpdate1Proc
#endif


#ifdef MPI
    subroutine testcaseHermitianmatrixIndexedUpdate4Proc
      Implicit None

      CALL set_test_name ('Trying to update whole matrix')
      CALL testHermitianmatrixIndexedUpdate4Proc

    end subroutine testcaseHermitianmatrixIndexedUpdate4Proc
#endif

!------------------------------------------------------------------------------
! test testNewComplexMatrixserial
!------------------------------------------------------------------------------
    subroutine testNewComplexMatrixSerial
      implicit none

      integer, parameter   :: nrows = 9, ncols = 10

      Type (ComplexMatrix) :: matrix
      
        Call newmatrix(matrix, nrows, ncols)

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

      Call deletematrix(matrix)

    end subroutine testNewComplexMatrixSerial

!------------------------------------------------------------------------------
! test testNewComplexMatrix1proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewComplexMatrix1proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter   :: nrows = 9, ncols = 10
      Type (ComplexMatrix) :: matrix
      
      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call newmatrix(matrix, nrows, ncols, DISTRIBUTE_2D)

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        Call deletematrix(matrix)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      endif

    end subroutine testNewComplexMatrix1proc
#endif


!------------------------------------------------------------------------------
! test testNewComplexMatrix4proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewComplexMatrix4proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter   :: nrows = 9, ncols = 10
      Type (ComplexMatrix) :: matrix

      integer              :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newmatrix(matrix, nrows, ncols, DISTRIBUTE_2D)
        
        Select Case (MPIglobal%myprocrow)
          Case (0)
            nrows_loc = 5
          Case (1)
            nrows_loc = 4
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            ncols_loc = 6
          Case (1)
            ncols_loc = 4
        End Select

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        Call deletematrix(matrix)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testNewComplexMatrix4proc
#endif


! !------------------------------------------------------------------------------
! ! test testNewComplexMatrix4procDistRows
! !------------------------------------------------------------------------------
! #ifdef MPI
!     subroutine testNewComplexMatrix4procDistRows
!       implicit none
! 
!       integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
! 
!       integer, parameter   :: nrows = 9, ncols = 10
!       Type (ComplexMatrix) :: matrix
! 
!       integer              :: nrows_loc, ncols_loc
! 
!       n_proc_rows_test = 4
!       n_proc_cols_test = 1
!       n_procs_test = n_proc_rows_test*n_proc_cols_test
!       call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
!       MPIglobal_1D%blocksize = 2
! 
!       if (MPIglobal_1D%rank < n_procs_test) then
!         Call getBlacsGridInfo(MPIglobal_1D)
! 
!         Call newmatrix(matrix, nrows, ncols, DISTRIBUTE_ROWS)
! 
!         ncols_loc = ncols
!         select case (MPIglobal%rank)
!           case (0)
!             nrows_loc = 3
!           case (1)
!             nrows_loc = 2
!           case (2)
!             nrows_loc = 2
!           case (3)
!             nrows_loc = 2
!         end select
! 
!         Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
!         Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
!         Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
!         Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
!         Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
!         Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
!         Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')
! 
!         Call deletematrix(matrix)
!         Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
!       end if
! 
!     end subroutine testNewComplexMatrix4procDistRows
! #endif


!------------------------------------------------------------------------------
! test testNewComplexMatrix4procDistCols
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewComplexMatrix4procDistCols
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter   :: nrows = 9, ncols = 10
      Type (ComplexMatrix) :: matrix

      integer              :: nrows_loc, ncols_loc

      n_proc_rows_test = 1
      n_proc_cols_test = 4
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

        Call newmatrix(matrix, nrows, ncols, DISTRIBUTE_COLS)
        
        nrows_loc = nrows
        select case (MPIglobal%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        Call deletematrix(matrix)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testNewComplexMatrix4procDistCols
#endif


!------------------------------------------------------------------------------
! test testNewHermitianMatrixserial
!------------------------------------------------------------------------------
    subroutine testNewHermitianMatrixSerial
      implicit none

      integer, parameter     :: nmatp = 10
      Type (HermitianMatrix) :: matrix
      
      Call newmatrix(matrix, nmatp)

      Call assert_equals(nmatp, matrix%size, 'checking newmatrix size')
      Call assert_equals(nmatp, size(matrix%za,1), 'checking newmatrix size rows')
      Call assert_equals(nmatp, size(matrix%za,2), 'checking newmatrix size cols')
      Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

      Call deletematrix(matrix)

    end subroutine testNewHermitianMatrixSerial

!------------------------------------------------------------------------------
! test testNewHermitianMatrix1proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewHermitianMatrix1proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter        :: nmatp = 10
      integer, dimension(nmatp) :: loc_idx
      integer                   :: i
      Type (HermitianMatrix)    :: matrix
      
      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        loc_idx = (/(i,i=1,nmatp)/)

        Call newmatrix(matrix, nmatp)

        Call assert_equals(nmatp, matrix%size, 'checking newmatrix size')
        Call assert_equals(nmatp, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(nmatp, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')
        Call assert_equals(loc_idx, matrix%my_rows_idx, nmatp, 'checking newmatrix local rows')
        Call assert_equals(loc_idx, matrix%my_cols_idx, nmatp, 'checking newmatrix local cols')

        Call deletematrix(matrix)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      endif

    end subroutine testNewHermitianMatrix1proc
#endif


!------------------------------------------------------------------------------
! test testNewHermitianMatrix4proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewHermitianMatrix4proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter     :: nmatp = 10
      Type (HermitianMatrix) :: matrix

      integer                            :: nrows_loc, ncols_loc
      integer, dimension(:), allocatable :: loc_rows_idx, loc_cols_idx

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Select Case (MPIglobal%myprocrow)
          Case (0)
            nrows_loc = 6
          Case (1)
            nrows_loc = 4
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            ncols_loc = 6
          Case (1)
            ncols_loc = 4
        End Select

        Allocate(loc_rows_idx(nrows_loc), loc_cols_idx(ncols_loc))
        Select Case (MPIglobal%myprocrow)
          Case (0)
            loc_rows_idx = (/1,2,5,6,9,10/)
          Case (1)
            loc_rows_idx = (/3,4,7,8/)
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            loc_cols_idx = (/1,2,5,6,9,10/)
          Case (1)
            loc_cols_idx = (/3,4,7,8/)
        End Select

        Call newmatrix(matrix, nmatp, DISTRIBUTE_2D)
        
        Call assert_equals(nmatp, matrix%size, 'checking newmatrix size')
        Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')
        Call assert_equals(loc_rows_idx, matrix%my_rows_idx, nrows_loc, 'checking newmatrix local rows')
        Call assert_equals(loc_cols_idx, matrix%my_cols_idx, ncols_loc, 'checking newmatrix local cols')

        Call deletematrix(matrix)
        Deallocate(loc_rows_idx, loc_cols_idx)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testNewHermitianMatrix4proc
#endif

! not needed, so not tested
! !------------------------------------------------------------------------------
! ! test testNewHermitianMatrix4procDistRows
! !------------------------------------------------------------------------------
! #ifdef MPI
!     subroutine testNewHermitianMatrix4procDistRows
!       implicit none
! 
!       integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
! 
!       integer, parameter     :: nmatp = 10
!       Type (HermitianMatrix) :: matrix
! 
!       integer                            :: nrows_loc, ncols_loc, i
!       integer, dimension(:), allocatable :: loc_rows_idx, loc_cols_idx
! 
!       n_proc_rows_test = 4
!       n_proc_cols_test = 1
!       n_procs_test = n_proc_rows_test*n_proc_cols_test
!       call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
!       MPIglobal_1D%blocksize = 2
! 
!       if (MPIglobal_1D%rank < n_procs_test) then
!         Call getBlacsGridInfo(MPIglobal_1D)
! 
!         ncols_loc = nmatp
!         select case (MPIglobal_1D%rank)
!           case (0)
!             nrows_loc = 4
!           case (1)
!             nrows_loc = 2
!           case (2)
!             nrows_loc = 2
!           case (3)
!             nrows_loc = 2
!         end select
! 
!         Allocate(loc_rows_idx(nrows_loc), loc_cols_idx(nmatp))
!         Select Case (MPIglobal_1D%myprocrow)
!           Case (0)
!             loc_rows_idx = (/1,2,9,10/)
!           Case (1)
!             loc_rows_idx = (/3,4/)
!           Case (2)
!             loc_rows_idx = (/5,6/)
!           Case (3)
!             loc_rows_idx = (/7,8/)
!         End Select
!         loc_cols_idx = (/(i,i=1,nmatp)/)
! 
!         Call newmatrix(matrix, nmatp, DISTRIBUTE_ROWS)
! 
!         Call assert_equals(nmatp, matrix%size, 'checking newmatrix size')
!         Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
!         Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
!         Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
!         Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
!         Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')
!         Call assert_equals(loc_rows_idx, matrix%my_rows_idx, nrows_loc, 'checking newmatrix local rows')
!         Call assert_equals(loc_cols_idx, matrix%my_cols_idx, ncols_loc, 'checking newmatrix local cols')
! 
!         Call deletematrix(matrix)
!         Deallocate(loc_rows_idx, loc_cols_idx)
!         Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
!       end if
! 
!     end subroutine testNewHermitianMatrix4procDistRows
! #endif


!------------------------------------------------------------------------------
! test testNewHermitianMatrix4procDistCols
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewHermitianMatrix4procDistCols
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter     :: nmatp = 10
      Type (HermitianMatrix) :: matrix

      integer                            :: nrows_loc, ncols_loc, i
      integer, dimension(:), allocatable :: loc_rows_idx, loc_cols_idx

      n_proc_rows_test = 1
      n_proc_cols_test = 4
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Allocate(loc_rows_idx(nmatp), loc_cols_idx(ncols_loc))
        loc_rows_idx = (/(i,i=1,nmatp)/)
        Select Case (MPIglobal_1D%myproccol)
          Case (0)
            loc_cols_idx = (/1,2,9,10/)
          Case (1)
            loc_cols_idx = (/3,4/)
          Case (2)
            loc_cols_idx = (/5,6/)
          Case (3)
            loc_cols_idx = (/7,8/)
        End Select

        Call newmatrix(matrix, nmatp, DISTRIBUTE_COLS)
        
        Call assert_equals(nmatp, matrix%size, 'checking newmatrix size')
        Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')
        Call assert_equals(loc_rows_idx, matrix%my_rows_idx, nrows_loc, 'checking newmatrix local rows')
        Call assert_equals(loc_cols_idx, matrix%my_cols_idx, ncols_loc, 'checking newmatrix local cols')

        Call deletematrix(matrix)
        Deallocate(loc_rows_idx, loc_cols_idx)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testNewHermitianMatrix4procDistCols
#endif


!------------------------------------------------------------------------------
! test testNewsystemserial
!------------------------------------------------------------------------------
    subroutine testNewSystemSerial
      implicit none

      integer, parameter :: nmatp = 9
      integer, parameter :: ngp   = 7
      Type (evsystem)    :: system
      
      Call newsystem (system, nmatp, ngp)

      Call assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
      Call assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
      Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking newsystem H size rows')
      Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking newsystem H size cols')
      Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')

      Call assert_equals(nmatp, system%overlap%size,  'checking newsystem S rank')
      Call assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
      Call assert_equals(nmatp, size(system%overlap%za,1), 'checking newsystem S size rows')
      Call assert_equals(nmatp, size(system%overlap%za,2), 'checking newsystem S size cols')
      Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')

      Call assert_equals(ngp, system%ngp, 'checking newsystem ngp')
      Call assert_equals(ngp, system%ngp_loc_rows, 'checking newsystem ngp_loc_rows')
      Call assert_equals(ngp, system%ngp_loc_cols, 'checking newsystem ngp_loc_cols')
      
      Call deletesystem(system)
    end subroutine testNewSystemSerial


!------------------------------------------------------------------------------
! test testNewsystem1proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewSystem1proc
      implicit none

      integer, parameter :: nmatp = 9
      integer, parameter :: ngp   = 7
      Type (evsystem)    :: system
      integer, dimension(nmatp) :: loc_idx
      
      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t, i

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        loc_idx = (/(i,i=1,nmatp)/)

        Call newsystem (system, nmatp, ngp)

        Call assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        Call assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        Call assert_equals(nmatp, size(system%hamilton%za,1), 'checking newsystem H size rows')
        Call assert_equals(nmatp, size(system%hamilton%za,2), 'checking newsystem H size cols')
        Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')
        Call assert_equals(loc_idx, system%hamilton%my_rows_idx, nmatp, 'checking newsystem H local rows')
        Call assert_equals(loc_idx, system%hamilton%my_cols_idx, nmatp, 'checking newsystem H local cols')

        Call assert_equals(nmatp, system%overlap%size,  'checking newsystem S rank')
        Call assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        Call assert_equals(nmatp, size(system%overlap%za,1), 'checking newsystem S size rows')
        Call assert_equals(nmatp, size(system%overlap%za,2), 'checking newsystem S size cols')
        Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')
        Call assert_equals(loc_idx, system%overlap%my_rows_idx, nmatp, 'checking newsystem H local rows')
        Call assert_equals(loc_idx, system%overlap%my_cols_idx, nmatp, 'checking newsystem H local cols')

        Call assert_equals(ngp, system%ngp, 'checking newsystem ngp')
        Call assert_equals(ngp, system%ngp_loc_rows, 'checking newsystem ngp_loc_rows')
        Call assert_equals(ngp, system%ngp_loc_cols, 'checking newsystem ngp_loc_cols')
      
        Call deletesystem(system)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      endif

    end subroutine testNewSystem1proc
#endif

!------------------------------------------------------------------------------
! test testNewsystem4proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewSystem4proc
      implicit none

      integer, parameter :: nmatp = 9
      integer, parameter :: ngp   = 7
      Type (evsystem)    :: system

      integer :: nrows_loc, ncols_loc, ngp_rows_loc, ngp_cols_loc
      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      integer, dimension(:), allocatable :: loc_rows_idx, loc_cols_idx

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newsystem (system, nmatp, ngp)
        
        Select Case (MPIglobal%myprocrow)
          Case (0)
            nrows_loc    = 5
            ngp_rows_loc = 4
          Case (1)
            nrows_loc    = 4
            ngp_rows_loc = 3
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            ncols_loc = 5
            ngp_cols_loc = 4
          Case (1)
            ncols_loc = 4
            ngp_cols_loc = 3
        End Select

        Allocate(loc_rows_idx(nrows_loc), loc_cols_idx(ncols_loc))
        Select Case (MPIglobal%myprocrow)
          Case (0)
            loc_rows_idx = (/1,2,5,6,9,10/)
          Case (1)
            loc_rows_idx = (/3,4,7,8/)
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            loc_cols_idx = (/1,2,5,6,9,10/)
          Case (1)
            loc_cols_idx = (/3,4,7,8/)
        End Select

        Call assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        Call assert_equals(nrows_loc, system%hamilton%nrows_loc, 'checking newsystem H nrows_loc')
        Call assert_equals(ncols_loc, system%hamilton%ncols_loc, 'checking newsystem H ncols_loc')
        Call assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking newsystem H size rows')
        Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking newsystem H size cols')
        Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')
        Call assert_equals(loc_rows_idx, system%hamilton%my_rows_idx, nrows_loc, 'checking newsystem H local rows')
        Call assert_equals(loc_cols_idx, system%hamilton%my_cols_idx, ncols_loc, 'checking newsystem H local cols')

        Call assert_equals(nmatp, system%overlap%size, 'checking newsystem S rank')
        Call assert_equals(nrows_loc, system%overlap%nrows_loc, 'checking newsystem S nrows_loc')
        Call assert_equals(ncols_loc, system%overlap%ncols_loc, 'checking newsystem S ncols_loc')
        Call assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking newsystem S size rows')
        Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking newsystem S size cols')
        Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')
        Call assert_equals(loc_rows_idx, system%overlap%my_rows_idx, nrows_loc, 'checking newsystem S local rows')
        Call assert_equals(loc_cols_idx, system%overlap%my_cols_idx, ncols_loc, 'checking newsystem S local cols')

        Call assert_equals(ngp, system%ngp, 'checking newsystem ngp')
        Call assert_equals(ngp_rows_loc, system%ngp_loc_rows, 'checking newsystem ngp_loc_rows')
        Call assert_equals(ngp_cols_loc, system%ngp_loc_cols, 'checking newsystem ngp_loc_cols')
      
        Call deletesystem(system)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testNewSystem4proc
#endif


!------------------------------------------------------------------------------
! test testNewsystem4procDistCols
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewSystem4procDistCols
      implicit none

      integer, parameter :: nmatp = 11
      integer, parameter :: ngp   = 9
      Type (evsystem)    :: system

      integer :: nrows_loc, ncols_loc, ngp_rows_loc, ngp_cols_loc, i
      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      integer, dimension(:), allocatable :: loc_rows_idx, loc_cols_idx

      n_proc_rows_test = 1
      n_proc_cols_test = 4
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

        nrows_loc = nmatp
        ngp_rows_loc = ngp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc    = 4
            ngp_cols_loc = 3
          case (1)
            ncols_loc    = 3
            ngp_cols_loc = 2
          case (2)
            ncols_loc    = 2
            ngp_cols_loc = 2
          case (3)
            ncols_loc    = 2
            ngp_cols_loc = 2
        end select

        Allocate(loc_rows_idx(nmatp), loc_cols_idx(ncols_loc))
        loc_rows_idx = (/(i,i=1,nmatp)/)
        Select Case (MPIglobal_1D%myproccol)
          Case (0)
            loc_cols_idx = (/1,2,9,10/)
          Case (1)
            loc_cols_idx = (/3,4,11/)
          Case (2)
            loc_cols_idx = (/5,6/)
          Case (3)
            loc_cols_idx = (/7,8/)
        End Select

        Call newsystem (system, nmatp, ngp, DISTRIBUTE_COLS)
        
        Call assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        Call assert_equals(nrows_loc, system%hamilton%nrows_loc, 'checking newsystem H nrows_loc')
        Call assert_equals(ncols_loc, system%hamilton%ncols_loc, 'checking newsystem H ncols_loc')
        Call assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking newsystem H size rows')
        Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking newsystem H size cols')
        Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')
        Call assert_equals(loc_rows_idx, system%hamilton%my_rows_idx, nrows_loc, 'checking newsystem H local rows')
        Call assert_equals(loc_cols_idx, system%hamilton%my_cols_idx, ncols_loc, 'checking newsystem H local cols')

        Call assert_equals(nmatp, system%overlap%size, 'checking newsystem S rank')
        Call assert_equals(nrows_loc, system%overlap%nrows_loc, 'checking newsystem S nrows_loc')
        Call assert_equals(ncols_loc, system%overlap%ncols_loc, 'checking newsystem S ncols_loc')
        Call assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking newsystem S size rows')
        Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking newsystem S size cols')
        Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')
        Call assert_equals(loc_rows_idx, system%overlap%my_rows_idx, nrows_loc, 'checking newsystem S local rows')
        Call assert_equals(loc_cols_idx, system%overlap%my_cols_idx, ncols_loc, 'checking newsystem S local cols')

        Call assert_equals(ngp, system%ngp, 'checking newsystem ngp')
        Call assert_equals(ngp_rows_loc, system%ngp_loc_rows, 'checking newsystem ngp_loc_rows')
        Call assert_equals(ngp_cols_loc, system%ngp_loc_cols, 'checking newsystem ngp_loc_cols')
      
        Call deletesystem(system)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testNewSystem4procDistCols
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrixSerial_AxI
! testing multiplication C=(A**H)xI   
! I is a rectangular matrix, upper square part filled with identity
! general A
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxI
      implicit none

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B

      Call newmatrix(C, nmatp)
      Call newMatrix(A, nrowsf, nmatp)
      Call fillComplex(A%za, nrowsf, nmatp)
      Call newMatrix(B, nrowsf, nmatp)
      Call fillIdentity(B%za, nmatp)

      Call newmatrix(ref, nmatp)
      ref%za= conjg(transpose(A%za(1:nmatp,:)))

      Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

      Call assert_equals(nmatp, C%size, 'checking result rank')
      Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
      Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
      Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

      Call deletematrix(A)
      Call deletematrix(B)
      Call deletematrix(C)
      Call deletematrix(ref)

    end subroutine testHermitianMatrixMatrixSerial_AxI

!------------------------------------------------------------------------------
! test testHermitianMatrixMatrixSerial_AxB
! testing multiplication C=AxB   
! hermitian A
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxB
      implicit none

      integer, parameter :: nmatp = 9, nrowsf = 10
      integer            :: row, col

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      
      Call newmatrix(C, nmatp)

      Call newMatrix(A, nrowsf, nmatp)
      Call fillComplex(A%za, nrowsf, nmatp)

      ! create matrix of type
      !        /  0  0  i  \
      !       |   0  i  0  |
      !       |   i  0  0  |
      !       |   0  0  0  |
      !        \  0  0  i /   -> one additional entry in last row
      Call newMatrix(B, nrowsf, nmatp)
      do row=1,nmatp
        B%za(row, nmatp-row+1) = cmplx(0,1,8)
      end do
      B%za(nrowsf,nmatp) = cmplx(0,1,8)

      Call newmatrix(ref, nmatp)
      do col=1,nmatp
        ref%za(:,col) = conjg(A%za(nmatp-col+1,:))*cmplx(0,1,8)
      end do
      ref%za(:,nmatp) = ref%za(:,nmatp) + conjg(A%za(nrowsf,:))*cmplx(0,1,8)

      Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

      Call assert_equals(nmatp, C%size, 'checking result rank')
      Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
      Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
      Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

      Call deletematrix(A)
      Call deletematrix(B)
      Call deletematrix(C)
      Call deletematrix(ref)
    end subroutine testHermitianMatrixMatrixSerial_AxB


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrixSerial_AxB_square
! testing multiplication C=AxB with 
! general matrices and A=B^dagger
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxB_square
      implicit none

      integer, parameter :: nmatp = 9
      integer            :: row, col
      double precision :: pi
      parameter (pi=3.1415926535897932385d0)
      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B

      Call newmatrix(C, nmatp)
      Call newMatrix(A, nmatp, nmatp)
      Call newMatrix(B, nmatp, nmatp)

      do row=1,nmatp
       do col=1,nmatp
        A%za(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
       enddo
      end do
      B%za=A%za

      Call newmatrix (ref, nmatp)
      ref%za(:,:)=cmplx(0d0,0d0,8)
      do row=1,nmatp
        ref%za(row,row) = cmplx(dble(nmatp),0d0,8)
      end do
      Call HermitianMatrixMatrix(C,A,B,nmatp,nmatp,nmatp)
!     write(*,*) A%za
!      stop
      CALL assert_equals(nmatp, C%size, 'checking result rank')
      CALL assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
      CALL assert_equals(nmatp, size(C%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(C%za,2), 'checking result size cols')
      CALL assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

      CALL deletematrix(A)
      CALL deletematrix(B)
      CALL deletematrix(C)
      CALL deletematrix(ref)
    end subroutine testHermitianMatrixMatrixSerial_AxB_square


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrixSerial_IxIpC
! testing multiplication C=IxI+C   
! I is a rectangular matrix, upper square part filled with identity
! complex C
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_IxIpC
      implicit none

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      integer :: i
      
      Call newmatrix(C, nmatp)
      Call fillHermitian(C%za, nmatp)
      Call newMatrix(A, nrowsf, nmatp)
      Call fillIdentity(A%za, nmatp)
      Call newMatrix(B, nrowsf, nmatp)
      Call fillIdentity(B%za, nmatp)

      Call newmatrix(ref, nmatp)
      ref%za = C%za
      do i=1,nmatp
        ref%za(i,i) = ref%za(i,i)+cmplx(1,0,8)
      end do

      Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

      Call assert_equals(nmatp, C%size, 'checking result rank')
      Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
      Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
      Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
      Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

      Call deletematrix(A)
      Call deletematrix(B)
      Call deletematrix(C)
      Call deletematrix(ref)
    end subroutine testHermitianMatrixMatrixSerial_IxIpC


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrixSerial_AxBpC_bighamiltonian
! testing multiplication C=AxB+C   
! C is bigger than AxB, so the multiplication should work on a submatrix of C
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxBpC_bighamiltonian
      implicit none

      integer, parameter :: nmatp = 10, ngp = 8, nrowsf = 4

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      integer                             :: i
      
      ! fill only the submatrix with the hermitian stuff, the rest set to 1
      Call newmatrix(C, nmatp)
      C%za = cmplx(1,0,8)
      Call fillHermitian(C%za, ngp)

      Call newMatrix(A, nrowsf, ngp)
      Call fillComplex(A%za, nrowsf, ngp)

      ! create matrix 
      !        /  1  0  0  0  i  0  0  0  \
      !        /  0  1  0  0  0  i  0  0  \
      !        /  0  0  1  0  0  0  i  0  \
      !        \  0  0  0  1  0  0  0  i /   
      Call newMatrix(B, nrowsf, ngp)
      do i=1,nrowsf
        B%za(i,i)   = cmplx(1,0,8)
        B%za(i,i+4) = cmplx(0,1,8)
      end do

      ! the product affects only the submatrix, the rest has to be unchanged
      Call newmatrix(ref, nmatp)
      ref%za(1:8,1:4) = conjg(transpose(A%za))
      ref%za(1:8,5:8) = conjg(transpose(A%za))*cmplx(0,1,8)
      ref%za          = ref%za + C%za

      Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,ngp)

      Call assert_equals(ref%za, C%za, ngp, ngp, tol, 'checking result numbers')

      Call deletematrix(A)
      Call deletematrix(B)
      Call deletematrix(C)
      Call deletematrix(ref)

    end subroutine testHermitianMatrixMatrixSerial_AxBpC_bighamiltonian


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxI
! testing multiplication C=(A**H)xI   in MPI mode with 1 proc
! I is a rectangular matrix, upper square part filled with identity
! general A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxI
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B


      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call newmatrix(C, nmatp)

        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
        Call fillComplex(A%za, nrowsf, nmatp)

        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
        Call fillIdentity(B%za, nmatp)

        Call newmatrix(ref, nmatp)
        ref%za= conjg(transpose(A%za(1:nmatp,:)))

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxI
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxB
! testing multiplication C=AxB   in MPI mode with 1 proc
! hermitian A, complex B
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxB
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10
      integer            :: row, col

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B


      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call newmatrix(C, nmatp)

        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
        Call fillComplex(A%za, nrowsf, nmatp)

        ! create matrix of type
        !        /  0  0  i  \
        !       |   0  i  0  |
        !       |   i  0  0  |
        !       |   0  0  0  |
        !        \  0  0  i /   -> one additional entry in last row
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
        do row=1,nmatp
          B%za(row, nmatp-row+1) = cmplx(0,1,8)
        end do
        B%za(nrowsf,nmatp) = cmplx(0,1,8)

        Call newmatrix(ref, nmatp)
        do col=1,nmatp
          ref%za(:,col) = conjg(A%za(nmatp-col+1,:))*cmplx(0,1,8)
        end do
        ref%za(:,nmatp) = ref%za(:,nmatp) + conjg(A%za(nrowsf,:))*cmplx(0,1,8)

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxB
#endif

!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxB_square
! testing multiplication C=AxB   in MPI mode with 1 proc
! general matrices and A=B^dagger
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxB_square
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      integer            :: row, col
      double precision :: pi
      parameter (pi=3.1415926535897932385d0)
      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B


      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newmatrix(C, nmatp)

        Call newMatrix(A, nmatp, nmatp, DISTRIBUTE_2D)
        do row=1,nmatp
         do col=1,nmatp
          A%za(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
         enddo
        end do

        Call newMatrix(B, nmatp, nmatp, DISTRIBUTE_2D)
        B%za=A%za

        Call newmatrix (ref, nmatp)
        ref%za(:,:)=cmplx(0d0,0d0,8)
        do row=1,nmatp
         ref%za(row,row) = cmplx(dble(nmatp),0d0,8)
        end do


        Call HermitianMatrixMatrix(C,A,B,nmatp,nmatp,nmatp)

        CALL assert_equals(nmatp, C%size, 'checking result rank')
        CALL assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        CALL assert_equals(nmatp, size(C%za,1), 'checking result size rows')
        CALL assert_equals(nmatp, size(C%za,2), 'checking result size cols')
        CALL assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

        CALL deletematrix(A)
        CALL deletematrix(B)
        CALL deletematrix(C)
        CALL deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxB_square
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_IxIpC
! testing multiplication C=AxI+C   
! hermitian A
! setting C to a constant, except the last col which is 0
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_IxIpC
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      integer                :: i

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call newmatrix(C, nmatp)
        Call fillHermitian(C%za, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
        Call fillIdentity(A%za, nmatp)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
        Call fillIdentity(B%za, nmatp)

        Call newmatrix(ref, nmatp)
        ref%za = C%za
        do i=1,nmatp
          ref%za(i,i) = ref%za(i,i)+cmplx(1,0,8)
        end do

        Call HermitianMatrixMatrix(C,A,B,nmatp,nmatp,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nmatp, size(C%za,1), 'checking result size rows')
        Call assert_equals(nmatp, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix1Proc_IxIpC
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxBpC_bighamiltonian
! testing multiplication C=AxB+C   in MPI mode with 1 proc
! C is bigger than AxB, so the multiplication should work on a submatrix of C
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxBpC_bighamiltonian
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 10, ngp = 8, nrowsf = 4

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      integer                             :: nrows_loc, ncols_loc
      integer                             :: i

      n_procs_test = 1
      call setupProcGrid(1, 1, MPIglobal%comm,    MPIglobal%context,    ierror_t)
      MPIglobal%blocksize    = 2
      MPIglobal_1D%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        ! fill only the submatrix with the hermitian stuff, the rest set to 1
        Call newmatrix(C, nmatp)
        C%za = cmplx(1,0,8)
        Call fillHermitian(C%za, ngp)

        Call newMatrix(A, nrowsf, ngp, DISTRIBUTE_2D)
        Call fillComplex(A%za, nrowsf, ngp)

        ! create matrix 
        !        /  1  0  0  0  i  0  0  0  \
        !        /  0  1  0  0  0  i  0  0  \
        !        /  0  0  1  0  0  0  i  0  \
        !        \  0  0  0  1  0  0  0  i /   
        Call newMatrix(B, nrowsf, ngp, DISTRIBUTE_2D)
        do i=1,nrowsf
          B%za(i,i)   = cmplx(1,0,8)
          B%za(i,i+4) = cmplx(0,1,8)
        end do

        ! the product affects only the submatrix, the rest has to be unchanged
        Call newmatrix(ref, nmatp)
        ref%za(1:8,1:4) = conjg(transpose(A%za))
        ref%za(1:8,5:8) = conjg(transpose(A%za))*cmplx(0,1,8)
        ref%za          = ref%za + C%za

        nrows_loc = nmatp
        ncols_loc = nmatp

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,ngp)

        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxBpC_bighamiltonian
#endif




!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxI
! testing multiplication C=AxI   (I=identity) in MPI mode with 4 procs
! hermitian A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxI
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: ref_global
      integer                             :: nrows_loc, ncols_loc


      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newmatrix(A, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ref_global = conjg(transpose(A_global(1:nmatp,:)))
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 3
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        Call assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        Call assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxI
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxB
! testing multiplication C=AxB   in MPI mode with 4 procs
! complex A, complex B
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxB
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10
      integer            :: row, col

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: ref_global
      integer                             :: nrows_loc, ncols_loc


      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        ! create matrix of type
        !        /  0  0  i  \
        !       |   0  i  0  |
        !       |   i  0  0  |
        !       |   0  0  0  |
        !        \  0  0  i /   -> one additional entry in last row
        do row=1,nmatp
          B_global(row, nmatp-row+1) = cmplx(0,1,8)
        end do
        B_global(nrowsf,nmatp) = cmplx(0,1,8)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        do col=1,nmatp
          ref_global(:,col) = conjg(A_global(nmatp-col+1,:))*cmplx(0,1,8)
        end do
        ref_global(:,nmatp) = ref_global(:,nmatp) + conjg(A_global(nrowsf,:))*cmplx(0,1,8)
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 3
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        Call assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        Call assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxB
#endif

!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxB_square
! testing multiplication C=AxB   in MPI mode with 4 procs
! general matrices and A=B^dagger
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxB_square
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 9
      integer            :: row, col
      double precision :: pi
      parameter (pi=3.1415926535897932385d0)

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global, ref_global
      integer                            :: nrows_loc, ncols_loc

      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

        do row=1,nmatp
         do col=1,nmatp
          A_global(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
         enddo
        end do
        Call newMatrix(A, nmatp, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        B_global=A_global
        Call newMatrix(B, nmatp, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ref_global(:,:)=cmplx(0d0,0d0,8)
        do row=1,nmatp
         ref_global(row,row) = cmplx(dble(nmatp),0d0,8)
        end do
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 3
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call HermitianMatrixMatrix(C,A,B,nmatp,nmatp,nmatp)

        CALL assert_equals(nmatp, C%size, 'checking result rank')
        CALL assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        CALL assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        CALL assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        CALL assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        CALL assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        CALL assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        CALL deletematrix(A)
        CALL deletematrix(B)
        CALL deletematrix(C)
        CALL deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxB_square
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_IxIpC
! testing multiplication C=IxI+C   (I=identity) in MPI mode with 4 procs
! hermitian C
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_IxIpC
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: C_global, ref_global
      integer                             :: nrows_loc, ncols_loc
      integer                             :: i

      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call fillHermitian(C_global, nmatp)
        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal_1D)

        Call fillIdentity(A_global, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ref_global = C_global
        do i=1,nmatp
          ref_global(i,i) = ref_global(i,i)+cmplx(1,0,8)
        end do
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 3
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        Call assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        Call assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_IxIpC
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxBpC
! testing multiplication C=AxB+C   in MPI mode with 4 procs
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxBpC
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 8, nrowsf = 4

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: C_global, ref_global
      integer                             :: nrows_loc, ncols_loc
      integer                             :: i

      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call fillHermitian(C_global, nmatp)
        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal_1D)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        ! create matrix 
        !        /  1  0  0  0  i  0  0  0  \
        !        /  0  1  0  0  0  i  0  0  \
        !        /  0  0  1  0  0  0  i  0  \
        !        \  0  0  0  1  0  0  0  i /   
        do i=1,nrowsf
          B_global(i,i)   = cmplx(1,0,8)
          B_global(i,i+4) = cmplx(0,1,8)
        end do
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ref_global(:,1:4) = conjg(transpose(A_global))
        ref_global(:,5:8) = conjg(transpose(A_global))*cmplx(0,1,8)
        ref_global        = ref_global + C_global
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        ncols_loc = 2

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        Call assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        Call assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxBpC_bighamiltonian
! testing multiplication C=AxB+C   in MPI mode with 4 procs
! C is bigger than AxB, so the multiplication should work on a submatrix of C
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxBpC_bighamiltonian
      implicit none

      integer :: n_procs_test, ierror_t

      integer, parameter :: nmatp = 10, ngp = 8, nrowsf = 4

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,ngp)   :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: C_global, ref_global
      integer                             :: nrows_loc, ncols_loc
      integer                             :: i

      n_procs_test = 4
      call setupProcGrid(1, n_procs_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)
      
        ! fill only the submatrix with the hermitian stuff, the rest set to 1
        C_global = cmplx(1,0,8)
        Call fillHermitian(C_global, ngp)
        Call newmatrix(C, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal_1D)

        Call fillComplex(A_global, nrowsf, ngp)
        Call newMatrix(A, nrowsf, ngp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal_1D)

        ! create matrix 
        !        /  1  0  0  0  i  0  0  0  \
        !        /  0  1  0  0  0  i  0  0  \
        !        /  0  0  1  0  0  0  i  0  \
        !        \  0  0  0  1  0  0  0  i /   
        do i=1,nrowsf
          B_global(i,i)   = cmplx(1,0,8)
          B_global(i,i+4) = cmplx(0,1,8)
        end do
        Call newMatrix(B, nrowsf, ngp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ! the product affects only the submatrix, the rest has to be unchanged
        ref_global(1:8,1:4) = conjg(transpose(A_global))
        ref_global(1:8,5:8) = conjg(transpose(A_global))*cmplx(0,1,8)
        ref_global          = ref_global + C_global
        Call newmatrix(ref, nmatp, DISTRIBUTE_COLS)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal_1D)

        nrows_loc = nmatp
        select case (MPIglobal_1D%rank)
          case (0)
            ncols_loc = 4
          case (1)
            ncols_loc = 2
          case (2)
            ncols_loc = 2
          case (3)
            ncols_loc = 2
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,ngp)

        Call assert_equals(nmatp, C%size, 'checking result rank')
        Call assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
        Call assert_equals(nrows_loc, C%nrows_loc, 'checking result nrows_loc')
        Call assert_equals(ncols_loc, C%ncols_loc, 'checking result ncols_loc')
        Call assert_equals(nrows_loc, size(C%za,1), 'checking result size rows')
        Call assert_equals(ncols_loc, size(C%za,2), 'checking result size cols')
        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC_bighamiltonian
#endif


!------------------------------------------------------------------------------
! test testHermitianmatrixIndexedUpdateSerial
!------------------------------------------------------------------------------
    subroutine testHermitianmatrixIndexedUpdateSerial
      implicit none

      integer, parameter   :: nmatp=5
      complex(8), parameter :: I = cmplx(1,1.1D0,8), O = cmplx(0,0,8), C = cmplx(2,3,8)

      Type (HermitianMatrix) :: matrix
      complex(8), dimension(nmatp,nmatp) :: matrix_ref
      integer :: row, col

      matrix_ref(1,:) = (/ I, I, I, I, I /)
      matrix_ref(2,:) = (/ O, I, I, I, I /)
      matrix_ref(3,:) = (/ O, O, I, I, I /)
      matrix_ref(4,:) = (/ O, O, O, I, I /)
      matrix_ref(5,:) = (/ O, O, O, O, I /)
      matrix_ref = matrix_ref + C
      
      Call newmatrix(matrix, nmatp)
      matrix%za = C

      do row=1,nmatp
        do col=1,nmatp
          call Hermitianmatrix_indexedupdate(matrix, col, row, I)
        end do
      end do

!      Call assert_equals(matrix_ref, matrix%za, nrows_loc, ncols_loc, tol, 'checking result numbers')
      Call assert_equals(matrix_ref, matrix%za, nmatp, nmatp, 'checking result numbers')

      Call deletematrix(matrix)

    end subroutine testHermitianmatrixIndexedUpdateSerial


!------------------------------------------------------------------------------
! test testHermitianmatrixIndexedUpdate1Proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianmatrixIndexedUpdate1Proc
      implicit none

      integer, parameter   :: nmatp=5
      complex(8), parameter :: I = cmplx(1,1.1D0,8), O = cmplx(0,0,8), C = cmplx(2,3,8)

      Type (HermitianMatrix) :: matrix
      complex(8), dimension(nmatp,nmatp) :: matrix_ref
      integer :: row, col

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        matrix_ref(1,:) = (/ I, I, I, I, I /)
        matrix_ref(2,:) = (/ O, I, I, I, I /)
        matrix_ref(3,:) = (/ O, O, I, I, I /)
        matrix_ref(4,:) = (/ O, O, O, I, I /)
        matrix_ref(5,:) = (/ O, O, O, O, I /)
        matrix_ref = matrix_ref + C
      
        Call newmatrix(matrix, nmatp)
        matrix%za = C

        do row=1,nmatp
          do col=1,nmatp
            call Hermitianmatrix_indexedupdate(matrix, col, row, I)
          end do
        end do

        Call assert_equals(matrix_ref, matrix%za, nmatp, nmatp, 'checking result numbers')

        Call deletematrix(matrix)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      endif

    end subroutine testHermitianmatrixIndexedUpdate1Proc
#endif


!------------------------------------------------------------------------------
! test testHermitianmatrixIndexedUpdate4Proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianmatrixIndexedUpdate4Proc
      implicit none

      integer, parameter    :: nmatp=5
      complex(8), parameter :: I = cmplx(1,1.1D0,8), O = cmplx(0,0,8), C = cmplx(2,3,8)

      Type (HermitianMatrix)             :: matrix_ref, matrix
      complex(8), dimension(nmatp,nmatp) :: matrix_ref_global
      integer :: row, col

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t
      integer :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      Call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        matrix_ref_global(1,:) = (/ I, I, I, I, I /)
        matrix_ref_global(2,:) = (/ O, I, I, I, I /)
        matrix_ref_global(3,:) = (/ O, O, I, I, I /)
        matrix_ref_global(4,:) = (/ O, O, O, I, I /)
        matrix_ref_global(5,:) = (/ O, O, O, O, I /)
        matrix_ref_global = matrix_ref_global + C
        Call newmatrix(matrix_ref, nmatp)
        Call getBlockDistributedLoc(matrix_ref_global, matrix_ref%za, MPIglobal)

        Select Case (MPIglobal%myprocrow)
          Case (0)
            nrows_loc = 3
          Case (1)
            nrows_loc = 2
        End Select
        Select Case (MPIglobal%myproccol)
          Case (0)
            ncols_loc = 3
          Case (1)
            ncols_loc = 2
        End Select

        Call newmatrix(matrix, nmatp)
        matrix%za = C
        do row=1,nrows_loc
          do col=1,ncols_loc
            call Hermitianmatrix_indexedupdate(matrix, col, row, I)
          end do
        end do

        Call assert_equals(matrix_ref%za, matrix%za, nrows_loc, ncols_loc, 'checking result numbers')

        Call deletematrix(matrix)
        Call deletematrix(matrix_ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      endif

    end subroutine testHermitianmatrixIndexedUpdate4Proc
#endif

end module modfvsystem_test
