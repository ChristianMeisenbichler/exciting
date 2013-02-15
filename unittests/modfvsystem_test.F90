module modfvsystem_test
    use fruit
    use modfvsystem
#ifdef MPI
    use modmpi
#endif

    implicit none

    real(8), parameter :: tol  = 1E-13
    real(8), parameter :: zero = 0.0D0

    contains

#ifdef MPI
    subroutine setupProcGrid(nprocrows_sub, nproccols_sub, comm, context, ierror)
      implicit none

      integer, intent(inout) :: nprocrows_sub, nproccols_sub
      integer, intent(out)   :: comm, context, ierror
      integer :: nprocs, i, world_group, procs_group

      nprocs = nprocrows_sub*nproccols_sub
      Call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror)
      Call MPI_GROUP_INCL(world_group, nprocs, (/(i,i=0,nprocs-1)/), procs_group, ierror)
      Call MPI_COMM_CREATE(MPI_COMM_WORLD, procs_group, comm, ierror);

      Call BLACS_GET(0, 0, context)
      Call BLACS_GRIDINIT(context,'r', nprocrows_sub, nproccols_sub)

    end subroutine setupProcGrid
#endif


#ifdef MPI
    subroutine getBlacsGridInfo(MPIdata)
      implicit none
      Type(MPIinfo), intent(inout) :: MPIdata

      Call BLACS_GRIDINFO(MPIdata%context, MPIdata%nprocrows, MPIdata%nproccols, MPIdata%myprocrow, MPIdata%myproccol)

    end subroutine getBlacsGridInfo
#endif


    subroutine fillIdentity(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row
      mat = 0
      do row=1,size
        mat(row, row) = cmplx(1.0D0, 0.0D0,8)
      end do

    end subroutine fillIdentity


    subroutine fillHermitian(mat, size)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: size

      integer :: row, col

      do row=1,size
        do col=1,row  
          if (row .ne. col) then 
            mat(row, col) = cmplx(row, col,8)
            mat(col, row) = cmplx(row,-col,8)
          else
            mat(col, row) = cmplx(row, 0.0D0,8)
          end if
        end do
      end do
    end subroutine fillHermitian


    subroutine fillComplex(mat, nrows, ncols)
      implicit none

      complex(8), dimension(:,:), intent(out) :: mat
      integer, intent(in)                     :: nrows, ncols

      integer :: row, col

      do row=1,nrows
        do col=1,ncols  
          mat(row, col) = cmplx(row,col+10,8)
        end do
      end do
    end subroutine fillComplex




    
    ! extract the local data for the specified process in a grid
    ! no communication!! the global matrix has to be present locally!!
#ifdef MPI
    subroutine getBlockDistributedLoc(global, local, MPIdata)
      implicit none

      Integer, External   :: NUMROC

      ! arguments
      complex(8), dimension(:,:), intent(in)  :: global
      complex(8), dimension(:,:), intent(out) :: local
      Type(MPIinfo), intent(in)               :: MPIdata

      ! local variables
      integer :: n_procs_x, n_procs_y, proc_x, proc_y, blocksize
      integer :: n_rows_glob, n_cols_glob, n_rows_loc, n_cols_loc, &
                 n_blocks_x_loc, i_block_x, glob_x_start, glob_x_end, loc_x_start, loc_x_end, &
                 n_blocks_y_loc, i_block_y, glob_y_start, glob_y_end, loc_y_start, loc_y_end

      n_procs_x = MPIdata%nproccols
      n_procs_y = MPIdata%nprocrows
      proc_x    = MPIdata%myproccol
      proc_y    = MPIdata%myprocrow
      blocksize = MPIdata%blocksize

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

      integer, parameter   :: nrows = 9, ncols = 10

      Type (ComplexMatrix) :: matrix
      
        Call newmatrix (matrix, nrows, ncols)

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
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newmatrix (matrix, nrows, ncols)

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        Call deletematrix(matrix)
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

        Call newmatrix (matrix, nrows, ncols)
        
        select case (MPIglobal%rank)
          case (0)
            nrows_loc = 5
            ncols_loc = 6
          case (1)
            nrows_loc = 5
            ncols_loc = 4
          case (2)
            nrows_loc = 4
            ncols_loc = 6
          case (3)
            nrows_loc = 4
            ncols_loc = 4
        end select

        Call assert_equals(nrows, matrix%nrows, 'checking newmatrix nrows')
        Call assert_equals(ncols, matrix%ncols, 'checking newmatrix ncols')
        Call assert_equals(nrows_loc, matrix%nrows_loc, 'checking newmatrix nrows_loc')
        Call assert_equals(ncols_loc, matrix%ncols_loc, 'checking newmatrix ncols_loc')
        Call assert_equals(nrows_loc, size(matrix%za,1), 'checking newmatrix size rows')
        Call assert_equals(ncols_loc, size(matrix%za,2), 'checking newmatrix size cols')
        Call assert_equals(zero, sum(abs(matrix%za)), tol, 'checking newmatrix =0')

        Call deletematrix(matrix)
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

      Call deletesystem(system)
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
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newsystem (system, nmatp)

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

        Call deletesystem(system)
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
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)

        Call newsystem (system, nmatp)
        
        select case (MPIglobal%rank)
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

        Call assert_equals(nmatp, system%hamilton%size, 'checking newsystem H rank')
        Call assert_equals(nrows_loc, system%hamilton%nrows_loc, 'checking newsystem H nrows_loc')
        Call assert_equals(ncols_loc, system%hamilton%ncols_loc, 'checking newsystem H ncols_loc')
        Call assert_true(.not. system%hamilton%ludecomposed, 'checking newsystem H ludecomposed')
        Call assert_equals(nrows_loc, size(system%hamilton%za,1), 'checking newsystem H size rows')
        Call assert_equals(ncols_loc, size(system%hamilton%za,2), 'checking newsystem H size cols')
        Call assert_equals(zero, sum(abs(system%hamilton%za)), tol, 'checking newsystem H=0')

        Call assert_equals(nmatp, system%overlap%size, 'checking newsystem S rank')
        Call assert_equals(nrows_loc, system%overlap%nrows_loc, 'checking newsystem S nrows_loc')
        Call assert_equals(ncols_loc, system%overlap%ncols_loc, 'checking newsystem S ncols_loc')
        Call assert_true(.not. system%overlap%ludecomposed, 'checking newsystem S ludecomposed')
        Call assert_equals(nrows_loc, size(system%overlap%za,1), 'checking newsystem S size rows')
        Call assert_equals(ncols_loc, size(system%overlap%za,2), 'checking newsystem S size cols')
        Call assert_equals(zero, sum(abs(system%overlap%za)), tol, 'checking newsystem S=0')

        Call deletesystem(system)
      end if

    end subroutine testNewSystem4proc
#endif


!------------------------------------------------------------------------------
! test testNewComplexMatrixserial_AxI
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

      Call newmatrix (ref, nmatp)
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
! test testNewComplexMatrixserial_AxB
! testing multiplication C=AxB   
! hermitian A
!------------------------------------------------------------------------------
    subroutine testHermitianMatrixMatrixSerial_AxB
      implicit none

      integer, parameter :: nmatp = 9, nrowsf = 10
      integer            :: row, col

      Type (HermitianMatrix) :: C, ref
      Type (ComplexMatrix)   :: A, B
      
      Call newmatrix (C, nmatp)

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
! test testNewComplexMatrixserial_AxB_square
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
! test testNewComplexMatrixserial_IxIpC
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

      Call newmatrix (ref, nmatp)
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
      
        Call newmatrix (C, nmatp)

        Call newMatrix(A, nrowsf, nmatp)
        Call fillComplex(A%za, nrowsf, nmatp)

        Call newMatrix(B, nrowsf, nmatp)
        Call fillIdentity(B%za, nmatp)

        Call newmatrix (ref, nmatp)
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
      
        Call newmatrix (C, nmatp)

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

        Call newmatrix (ref, nmatp)
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

        Call newmatrix (C, nmatp)

        Call newMatrix(A, nmatp, nmatp)
        do row=1,nmatp
         do col=1,nmatp
          A%za(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
         enddo
        end do

        Call newMatrix(B, nmatp, nmatp)
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
        Call newMatrix(A, nrowsf, nmatp)
        Call fillIdentity(A%za, nmatp)
        Call newMatrix(B, nrowsf, nmatp)
        Call fillIdentity(B%za, nmatp)

        Call newmatrix (ref, nmatp)
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
      end if

    end subroutine testHermitianMatrixMatrix1Proc_IxIpC
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

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: ref_global
      integer                             :: nrows_loc, ncols_loc


      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call newmatrix(C, nmatp)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newmatrix(A, nrowsf, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nrowsf, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        ref_global = conjg(transpose(A_global(1:nmatp,:)))
        Call newmatrix (ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
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

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10
      integer            :: row, col

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: ref_global
      integer                             :: nrows_loc, ncols_loc


      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call newmatrix (C, nmatp)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newMatrix(A, nrowsf, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

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
        Call newMatrix(B, nrowsf, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        do col=1,nmatp
          ref_global(:,col) = conjg(A_global(nmatp-col+1,:))*cmplx(0,1,8)
        end do
        ref_global(:,nmatp) = ref_global(:,nmatp) + conjg(A_global(nrowsf,:))*cmplx(0,1,8)
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
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

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      integer            :: row, col
      double precision :: pi
      parameter (pi=3.1415926535897932385d0)

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global, ref_global
      integer                            :: nrows_loc, ncols_loc


      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        do row=1,nmatp
         do col=1,nmatp
          A_global(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
         enddo
        end do
        B_global=A_global

        ref_global(:,:)=cmplx(0d0,0d0,8)
        do row=1,nmatp
         ref_global(row,row) = cmplx(dble(nmatp),0d0,8)
        end do

        Call getBlacsGridInfo(MPIglobal)

        Call newmatrix(C, nmatp)

        Call newMatrix(A, nmatp, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        Call newMatrix(B, nmatp, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        Call newmatrix (ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
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

    end subroutine testHermitianMatrixMatrix4Proc_AxB_square
#endif

!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_IxIpC
! testing multiplication C=AxI   (I=identity) in MPI mode with 4 procs
! hermitian A
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_IxIpC
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9
      complex, parameter :: addconst = (-1,1)

      Type (HermitianMatrix)             :: C, ref
      Type (ComplexMatrix)               :: A, B
      complex(8), dimension(nmatp,nmatp) :: A_global, B_global, C_global, ref_global
      integer                            :: nrows_loc, ncols_loc
      integer                            :: i

      n_proc_rows_test = 2
      n_proc_cols_test = 2
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal%comm, MPIglobal%context, ierror_t)
      MPIglobal%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        Call fillHermitian(C_global, nmatp)
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal)

        Call fillIdentity(A_global, nmatp)
        Call newMatrix(A, nmatp, nmatp)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nmatp, nmatp)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        ref_global = C_global
        do i=1,nmatp
          ref_global(i,i) = ref_global(i,i)+cmplx(1,0,8)
        end do
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
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
      end if

    end subroutine testHermitianMatrixMatrix4Proc_IxIpC
#endif




end module modfvsystem_test
