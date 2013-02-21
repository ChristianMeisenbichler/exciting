module modfvsystem_test
    use fruit
    use test_helpers
    use modfvsystem
#ifdef MPI
    use modmpi
#endif

    implicit none

    contains


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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testNewComplexMatrix4proc
#endif


!------------------------------------------------------------------------------
! test testNewComplexMatrix4procDistRows
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testNewComplexMatrix4procDistRows
      implicit none

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter   :: nrows = 9, ncols = 10
      Type (ComplexMatrix) :: matrix

      integer              :: nrows_loc, ncols_loc

      n_proc_rows_test = 4
      n_proc_cols_test = 1
      n_procs_test = n_proc_rows_test*n_proc_cols_test
      call setupProcGrid(n_proc_rows_test, n_proc_cols_test, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal_1D%blocksize = 2

      if (MPIglobal_1D%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal_1D)

        Call newmatrix(matrix, nrows, ncols, DISTRIBUTE_ROWS)

        ncols_loc = ncols
        select case (MPIglobal%rank)
          case (0)
            nrows_loc = 3
          case (1)
            nrows_loc = 2
          case (2)
            nrows_loc = 2
          case (3)
            nrows_loc = 2
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

    end subroutine testNewComplexMatrix4procDistRows
#endif


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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
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
        Call newmatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        ref_global = conjg(transpose(A_global(1:nmatp,:)))
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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
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
      
        Call newmatrix(C, nmatp)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
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
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
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
        Call getBlacsGridInfo(MPIglobal)

        do row=1,nmatp
         do col=1,nmatp
          A_global(row, col) = cmplx(cos(-dble(2*row*(col-1))/dble(nmatp)*pi),sin(-dble(2*row*(col-1))/dble(nmatp)*pi),8)
         enddo
        end do
        Call newMatrix(A, nmatp, nmatp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        B_global=A_global
        Call newMatrix(B, nmatp, nmatp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        ref_global(:,:)=cmplx(0d0,0d0,8)
        do row=1,nmatp
         ref_global(row,row) = cmplx(dble(nmatp),0d0,8)
        end do
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        Call newmatrix(C, nmatp)

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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
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

      integer :: n_procs_test, n_proc_rows_test, n_proc_cols_test, ierror_t

      integer, parameter :: nmatp = 9, nrowsf = 10

      Type (HermitianMatrix)              :: C, ref
      Type (ComplexMatrix)                :: A, B
      complex(8), dimension(nrowsf,nmatp) :: A_global, B_global
      complex(8), dimension(nmatp,nmatp)  :: C_global, ref_global
      integer                             :: nrows_loc, ncols_loc
      integer                             :: i

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
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        Call fillIdentity(B_global, nmatp)
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_2D)
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
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_IxIpC
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxBpC_coldist
! testing multiplication C=AxB+C   in MPI mode with 4 procs
! A and B are col-distributed, C in both dimensions
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxBpC_coldist
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
      call setupProcGrid(2, 2, MPIglobal%comm,    MPIglobal%context,    ierror_t)
      call setupProcGrid(1, 4, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize    = 2
      MPIglobal_1D%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call fillHermitian(C_global, nmatp)
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal)

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
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        nrows_loc = 4
        ncols_loc = 4

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC_coldist
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxBpC_rowdist
! testing multiplication C=AxB+C   in MPI mode with 4 procs
! A and B are row-distributed, C in both dimensions
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxBpC_rowdist
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
      call setupProcGrid(2, 2, MPIglobal%comm,    MPIglobal%context,    ierror_t)
      call setupProcGrid(4, 1, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize    = 2
      MPIglobal_1D%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)
      
        Call fillHermitian(C_global, nmatp)
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal)

        Call fillComplex(A_global, nrowsf, nmatp)
        Call newMatrix(A, nrowsf, nmatp, DISTRIBUTE_ROWS)
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
        Call newMatrix(B, nrowsf, nmatp, DISTRIBUTE_ROWS)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal_1D)

        ref_global(:,1:4) = conjg(transpose(A_global))
        ref_global(:,5:8) = conjg(transpose(A_global))*cmplx(0,1,8)
        ref_global        = ref_global + C_global
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        nrows_loc = 4
        ncols_loc = 4

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,nmatp)

        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC_rowdist
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
      call setupProcGrid(2, 2, MPIglobal%comm,    MPIglobal%context,    ierror_t)
      MPIglobal%blocksize    = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
      
        ! fill only the submatrix with the hermitian stuff, the rest set to 1
        C_global = cmplx(1,0,8)
        Call fillHermitian(C_global, ngp)
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal)

        Call fillComplex(A_global, nrowsf, ngp)
        Call newMatrix(A, nrowsf, ngp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(A_global, A%za, MPIglobal)

        ! create matrix 
        !        /  1  0  0  0  i  0  0  0  \
        !        /  0  1  0  0  0  i  0  0  \
        !        /  0  0  1  0  0  0  i  0  \
        !        \  0  0  0  1  0  0  0  i /   
        do i=1,nrowsf
          B_global(i,i)   = cmplx(1,0,8)
          B_global(i,i+4) = cmplx(0,1,8)
        end do
        Call newMatrix(B, nrowsf, ngp, DISTRIBUTE_2D)
        Call getBlockDistributedLoc(B_global, B%za, MPIglobal)

        ! the product affects only the submatrix, the rest has to be unchanged
        ref_global(1:8,1:4) = conjg(transpose(A_global))
        ref_global(1:8,5:8) = conjg(transpose(A_global))*cmplx(0,1,8)
        ref_global          = ref_global + C_global
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
          case (0)
            nrows_loc = 6
            ncols_loc = 6
          case (1)
            nrows_loc = 6
            ncols_loc = 4
          case (2)
            nrows_loc = 4
            ncols_loc = 6
          case (3)
            nrows_loc = 4
            ncols_loc = 4
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,ngp)

        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC_bighamiltonian
#endif


!------------------------------------------------------------------------------
! test testHermitianMatrixMatrix4Proc_AxBpC_coldist_bighamiltonian
! testing multiplication C=AxB+C   in MPI mode with 4 procs
! A and B are col-distributed, C in both dimensions
! C is bigger than AxB, so the multiplication should work on a submatrix of C
!------------------------------------------------------------------------------
#ifdef MPI
    subroutine testHermitianMatrixMatrix4Proc_AxBpC_coldist_bighamiltonian
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
      call setupProcGrid(2, 2, MPIglobal%comm,    MPIglobal%context,    ierror_t)
      call setupProcGrid(1, 4, MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      MPIglobal%blocksize    = 2
      MPIglobal_1D%blocksize = 2

      if (MPIglobal%rank < n_procs_test) then
        Call getBlacsGridInfo(MPIglobal)
        Call getBlacsGridInfo(MPIglobal_1D)
      
        ! fill only the submatrix with the hermitian stuff, the rest set to 1
        C_global = cmplx(1,0,8)
        Call fillHermitian(C_global, ngp)
        Call newmatrix(C, nmatp)
        Call getBlockDistributedLoc(C_global, C%za, MPIglobal)

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
        Call newmatrix(ref, nmatp)
        Call getBlockDistributedLoc(ref_global, ref%za, MPIglobal)

        select case (MPIglobal%rank)
          case (0)
            nrows_loc = 6
            ncols_loc = 6
          case (1)
            nrows_loc = 6
            ncols_loc = 4
          case (2)
            nrows_loc = 4
            ncols_loc = 6
          case (3)
            nrows_loc = 4
            ncols_loc = 4
        end select

        Call HermitianMatrixMatrix(C,A,B,nrowsf,nrowsf,ngp)

        Call assert_equals(ref%za, C%za, nrows_loc, ncols_loc, tol, 'checking result numbers')

        Call deletematrix(A)
        Call deletematrix(B)
        Call deletematrix(C)
        Call deletematrix(ref)
        Call finalizeProcGrid(MPIglobal%comm, MPIglobal%context, ierror_t)
        Call finalizeProcGrid(MPIglobal_1D%comm, MPIglobal_1D%context, ierror_t)
      end if

    end subroutine testHermitianMatrixMatrix4Proc_AxBpC_coldist_bighamiltonian
#endif


end module modfvsystem_test
