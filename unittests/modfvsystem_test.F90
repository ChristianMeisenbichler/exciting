module modfvsystem_test
    use fruit
    use modfvsystem
#ifdef MPI
    use modmpi
#endif

    implicit none

    real(8), parameter :: tol  = 1E-9
    real(8), parameter :: zero = 0.0D0

    contains

#ifdef MPI
    subroutine setupProcGrid(nprocrows_in, nproccols_in, comm, ierror)
      implicit none

      integer, intent(inout) :: nprocrows_in, nproccols_in
      integer, intent(out)   :: comm, ierror
      integer :: nprocs, i, world_group, procs_group

      nprocrows = nprocrows_in
      nproccols = nproccols_in
      nprocs = nprocrows*nproccols
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror)
      CALL MPI_GROUP_INCL(world_group, nprocs, (/(i,i=0,nprocs-1)/), procs_group, ierror)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD, procs_group, comm, ierror);

      CALL BLACS_GET(0, 0, context)
      CALL BLACS_GRIDINIT(context,'r', nprocrows, nproccols)

    end subroutine setupProcGrid
#endif


    subroutine setHermitian(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row, col

      do row=1,size
        do col=1,row  
          mat(row, col) = cmplx(row,col)
          if (row .ne. col) then 
            mat(col, row) = cmplx(row, -col)
          else
            mat(col, row) = cmplx(row, 0.0D0)
          end if
        end do
      end do
    end subroutine setHermitian


    subroutine setIdentity(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row
      mat = 0
      do row=1,size
        mat(row, row) = 1
      end do

    end subroutine setIdentity

    
    ! extract the local data for the specified process in a grid
    ! no communication!! the global matrix has to be present locally!!
#ifdef MPI
    subroutine getBlockDistributedLoc(global, local, n_procs_x, n_procs_y, proc_x, proc_y, blocksize)
      implicit none

      Integer, External   :: NUMROC

      ! arguments
      complex(8), dimension(:,:), intent(in)  :: global
      complex(8), dimension(:,:), intent(out) :: local
      integer,                    intent(in)  :: n_procs_x, n_procs_y, proc_x, proc_y, blocksize

      ! local variables
      integer :: n_rows_glob, n_cols_glob, n_rows_loc, n_cols_loc, &
                  n_blocks_x_loc, i_block_x, glob_x_start, glob_x_end, loc_x_start, loc_x_end, &
                  n_blocks_y_loc, i_block_y, glob_y_start, glob_y_end, loc_y_start, loc_y_end

      n_rows_glob = SIZE(global, 1)
      n_cols_glob = SIZE(global, 2)

      n_rows_loc = NUMROC(n_rows_glob, blocksize, proc_y, 0, n_procs_y)
      n_cols_loc = NUMROC(n_cols_glob, blocksize, proc_x, 0, n_procs_x)

      n_blocks_x_loc = CEILING(FLOAT(n_cols_loc)/blocksize)
      n_blocks_y_loc = CEILING(FLOAT(n_rows_loc)/blocksize)

      do i_block_x=0,n_blocks_x_loc-1
          glob_x_start = (proc_x + i_block_x*n_procs_x)*blocksize + 1
          glob_x_end   = MIN(glob_x_start + blocksize - 1, n_cols_glob)
          loc_x_start  = i_block_x*blocksize + 1
          loc_x_end    = MIN(loc_x_start + blocksize - 1, n_cols_loc)

          do i_block_y=0,n_blocks_y_loc-1
              glob_y_start = (proc_y + i_block_y*n_procs_y)*blocksize + 1
              glob_y_end   = MIN(glob_y_start + blocksize - 1, n_rows_glob)
              loc_y_start  = i_block_y*blocksize + 1
              loc_y_end    = MIN(loc_y_start + blocksize - 1, n_rows_loc)
              local(loc_y_start:loc_y_end,loc_x_start:loc_x_end) = global(glob_y_start:glob_y_end,glob_x_start:glob_x_end)
          end do
      end do

    end subroutine getBlockDistributedLoc
#endif



!------------------------------------------------------------------------------
! test testNewComplexMatrixserial
!------------------------------------------------------------------------------
    subroutine testNewComplexMatrixSerial
      implicit none

      integer, parameter :: nmatp = 9

      Type (ComplexMatrix) :: matrix
      
      Call newmatrix (matrix, nmatp)

      CALL assert_equals(nmatp, matrix%size, 'checking newmatrix rank')
      CALL assert_equals(nmatp, size(matrix%za,1), 'checking newmatrix size rows')
      CALL assert_equals(nmatp, size(matrix%za,2), 'checking newmatrix size cols')
      CALL assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

      CALL deletematrix(matrix)

    end subroutine testNewComplexMatrixSerial


!------------------------------------------------------------------------------
! test testNewComplexMatrix1proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewComplexMatrix1proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter   :: nmatp = 9
      Type (ComplexMatrix) :: matrix
      
      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

        Call newmatrix (matrix, nmatp)

        CALL assert_equals(nmatp, matrix%size, 'checking newmatrix H rank')
        CALL assert_equals(nmatp, size(matrix%za,1), 'checking newmatrix H size rows')
        CALL assert_equals(nmatp, size(matrix%za,2), 'checking newmatrix H size cols')
        CALL assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix H=0')

        CALL deletematrix(matrix)
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

      integer, parameter   :: nmatp = 9
      Type (ComplexMatrix) :: matrix

      integer              :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

        CALL newmatrix (matrix, nmatp)
        
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

        CALL assert_equals(nmatp, matrix%size, 'checking newmatrix rank')
        CALL assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        CALL assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        CALL assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        CALL assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        CALL assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        CALL deletematrix(matrix)
      end if

    end subroutine testNewComplexMatrix4proc
#endif


!------------------------------------------------------------------------------
! test testNewsystemserial
!------------------------------------------------------------------------------
    subroutine testNewSystemSerial
      implicit none

      integer, parameter :: nmatp = 9
      Type (evsystem) :: system
      
      Call newsystem (system, nmatp)

      CALL assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
      CALL assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
      CALL assert_equals(nmatp, size(system%hamilton%za,1), 'checking newsystem H size rows')
      CALL assert_equals(nmatp, size(system%hamilton%za,2), 'checking newsystem H size cols')
      CALL assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')

      CALL assert_equals(nmatp, system%overlap%size,  'checking newsystem S rank')
      CALL assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
      CALL assert_equals(nmatp, size(system%overlap%za,1), 'checking newsystem S size rows')
      CALL assert_equals(nmatp, size(system%overlap%za,2), 'checking newsystem S size cols')
      CALL assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')

      CALL deletesystem(system)
    end subroutine testNewSystemSerial

!------------------------------------------------------------------------------
! test testNewsystem1proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewSystem1proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      Type (evsystem) :: system
      
      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

        Call newsystem (system, nmatp)

        CALL assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        CALL assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        CALL assert_equals(nmatp, size(system%hamilton%za,1), 'checking newsystem H size rows')
        CALL assert_equals(nmatp, size(system%hamilton%za,2), 'checking newsystem H size cols')
        CALL assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')

        CALL assert_equals(nmatp, system%overlap%size,  'checking newsystem S rank')
        CALL assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        CALL assert_equals(nmatp, size(system%overlap%za,1), 'checking newsystem S size rows')
        CALL assert_equals(nmatp, size(system%overlap%za,2), 'checking newsystem S size cols')
        CALL assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')

        CALL deletesystem(system)
      endif

    end subroutine testNewSystem1proc
#endif

!------------------------------------------------------------------------------
! test testNewsystem4proc
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewSystem4proc
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      Type (evsystem)    :: system

      integer            :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nprocrows, nproccols, myprocrow, myproccol)

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

        CALL assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        CALL assert_equals(nrows_loc, system%hamilton%nrows_loc, 'checking newsystem H nrows_loc')
        CALL assert_equals(ncols_loc, system%hamilton%ncols_loc, 'checking newsystem H ncols_loc')
        CALL assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        CALL assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking newsystem H size rows')
        CALL assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking newsystem H size cols')
        CALL assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')

        CALL assert_equals(nmatp, system%overlap%size, 'checking newsystem S rank')
        CALL assert_equals(nrows_loc, system%overlap%nrows_loc, 'checking newsystem S nrows_loc')
        CALL assert_equals(ncols_loc, system%overlap%ncols_loc, 'checking newsystem S ncols_loc')
        CALL assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        CALL assert_equals(nrows_loc, size(system%overlap%za,1), 'checking newsystem S size rows')
        CALL assert_equals(ncols_loc, size(system%overlap%za,2), 'checking newsystem S size cols')
        CALL assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')

        CALL deletesystem(system)
      end if

    end subroutine testNewSystem4proc
#endif


!------------------------------------------------------------------------------
! test testNewComplexMatrixserial_AxI
! testing multiplication C=AxI   (I=identity)
! hermitian A
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxI
      implicit none

      integer, parameter :: nmatp = 9

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      
      Call newmatrix (C, nmatp)
      Call newMatrix(A, nmatp)
      Call setHermitian(A%za, nmatp)
      Call newMatrix(B, nmatp)
      Call setIdentity(B%za, nmatp)

      Call newmatrix (ref, nmatp)
      ref%za = A%za

      Call HermitianMatrixMatrixNew(C,A,B,nmatp,nmatp,nmatp)

      CALL assert_equals(nmatp, C%size, 'checking result rank')
      CALL assert_true(.not. C%ludecomposed, 'checking result ludecomposed')
      CALL assert_equals(nmatp, size(C%za,1), 'checking result size rows')
      CALL assert_equals(nmatp, size(C%za,2), 'checking result size cols')
      CALL assert_equals(ref%za, C%za, nmatp, nmatp, tol, 'checking result numbers')

      CALL deletematrix(A)
      CALL deletematrix(B)
      CALL deletematrix(C)
      CALL deletematrix(ref)
    end subroutine testHermitianMatrixMatrixSerial_AxI


!------------------------------------------------------------------------------
! test testNewComplexMatrixserial_AxB
! testing multiplication C=AxB   
! hermitian A
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxB
      implicit none

      integer, parameter :: nmatp = 9
      integer            :: row, col

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      
      Call newmatrix (C, nmatp)
      Call newMatrix(A, nmatp)
      Call setHermitian(A%za, nmatp)
      Call newMatrix(B, nmatp)
      do row=1,nmatp
        B%za(row, nmatp-row+1) = cmplx(0,1)
      end do

      Call newmatrix (ref, nmatp)
      do col=1,nmatp
        ref%za(:,col) = A%za(:,nmatp-col+1)*cmplx(0,1)
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
    end subroutine testHermitianMatrixMatrixSerial_AxB


!------------------------------------------------------------------------------
! test testNewComplexMatrixserial_AxIpc
! testing multiplication C=AxI+C   
! hermitian A
! setting C to a constant, except the last col which is 0
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxIpC
      implicit none

      integer, parameter :: nmatp = 9
      complex, parameter :: addconst = (-1,1)

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      
      Call newmatrix (C, nmatp)
      C%za(:,1:nmatp-1) = addconst
      Call newMatrix(A, nmatp)
      Call setHermitian(A%za, nmatp)
      Call newMatrix(B, nmatp)
      Call setIdentity(B%za, nmatp)

      Call newmatrix (ref, nmatp)
      ref%za(:,1:nmatp-1) = A%za(:,1:nmatp-1)+addconst
      ref%za(:,nmatp)     = A%za(:,nmatp)

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
    end subroutine testHermitianMatrixMatrixSerial_AxIpC


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxI
! testing multiplication C=AxI   (I=identity) in MPI mode with 1 proc
! hermitian A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxI
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B


      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
      
        Call newmatrix (C, nmatp)
        Call newMatrix(A, nmatp)
        Call setHermitian(A%za, nmatp)
        Call newMatrix(B, nmatp)
        Call setIdentity(B%za, nmatp)

        Call newmatrix (ref, nmatp)
        ref%za = A%za

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

      integer, parameter :: nmatp = 9
      integer            :: row, col

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B


      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
      
        Call newmatrix (C, nmatp)
        Call newMatrix(A, nmatp)
        Call setHermitian(A%za, nmatp)
        Call newMatrix(B, nmatp)
        do row=1,nmatp
          B%za(row, nmatp-row+1) = cmplx(0,1)
        end do

        Call newmatrix (ref, nmatp)
        do col=1,nmatp
          ref%za(:,col) = A%za(:,nmatp-col+1)*cmplx(0,1)
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
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxB
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxIpc
! testing multiplication C=AxI+C   
! hermitian A
! setting C to a constant, except the last col which is 0
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix1Proc_AxIpC
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      complex, parameter :: addconst = (-1,1)

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B

      n_proc_rows_test = 1
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
      
        Call newmatrix(C, nmatp)
        C%za(:,1:nmatp-1) = addconst
        Call newMatrix(A, nmatp)
        Call setHermitian(A%za, nmatp)
        Call newMatrix(B, nmatp)
        Call setIdentity(B%za, nmatp)

        Call newmatrix (ref, nmatp)
        ref%za(:,1:nmatp-1) = A%za(:,1:nmatp-1)+addconst
        ref%za(:,nmatp)     = A%za(:,nmatp)

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
      end if

    end subroutine testHermitianMatrixMatrix1Proc_AxIpC
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxI
! testing multiplication C=AxI   (I=identity) in MPI mode with 4 procs
! hermitian A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxI
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global
      integer                            :: nrows_loc, ncols_loc


      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nproccols, nproccols, myprocrow, myproccol)
      
        Call newmatrix (C, nmatp)

        Call setHermitian(A_global, nmatp)
        Call newMatrix(A, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        Call setIdentity(B_global, nmatp)
        Call newMatrix(B, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        Call newmatrix (ref, nmatp)
        ref%za = A%za

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
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxI
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix1Proc_AxI
! testing multiplication C=AxB   in MPI mode with 4 procs
! hermitian A, complex B
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxB
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      integer            :: row, col

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global, ref_global
      integer                            :: nrows_loc, ncols_loc


      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nproccols, nproccols, myprocrow, myproccol)
      
        Call newmatrix (C, nmatp)

        Call setHermitian(A_global, nmatp)
        Call newMatrix(A, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        do row=1,nmatp
          B_global(row, nmatp-row+1) = cmplx(0,1)
        end do
        Call newMatrix(B, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        do col=1,nmatp
          ref_global(:,col) = A_global(:,nmatp-col+1)*cmplx(0,1)
        end do
        Call newmatrix (ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

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
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxB
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxIpC
! testing multiplication C=AxI   (I=identity) in MPI mode with 4 procs
! hermitian A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxIpC
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      complex, parameter :: addconst = (-1,1)

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global, C_global, ref_global
      integer                            :: nrows_loc, ncols_loc

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, comm, ierror_t)
      blocksize = 2

      if (rank < n_procs_test) then
        CALL BLACS_GRIDINFO(context, nproccols, nproccols, myprocrow, myproccol)
      
        C_global(:,1:nmatp-1) = addconst
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        Call setHermitian(A_global, nmatp)
        Call newMatrix(A, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        Call setIdentity(B_global, nmatp)
        Call newMatrix(B, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

        ref_global(:,1:nmatp-1) = A_global(:,1:nmatp-1)+addconst
        ref_global(:,nmatp)     = A_global(:,nmatp)
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, nproccols, nprocrows, myproccol, myprocrow, blocksize)

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
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxIpC
#endif




end module modfvsystem_test
