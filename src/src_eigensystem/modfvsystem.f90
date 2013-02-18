! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE:  modfvsystem
! !DESCRIPTION:
! Module for setting up the eigensystem
! it is designed in a way that all other subroutines
! dealing with setting up and solving the system can acsess the
! data transparently allowing to choose from different datatypes
! more easily
!
!
! !REVISION HISTORY:
!   Created ??
!   Parallelized, Feb 2013 (G. Huhs - BSC)
!
!
!
  Module modfvsystem
! !USES:
#ifdef MPI
    use modmpi
#endif
    Implicit None

    !
    Type ComplexMatrix
       !
       Integer :: nrows, ncols 
       Integer :: nrows_loc, ncols_loc
       Complex (8), Pointer :: za (:, :)
#ifdef MPI
       Integer, Dimension(9)    :: desc
#endif
    End Type ComplexMatrix

    Type HermitianMatrix
       !
       Integer :: size   !TODO: in the end we should be able to get rid of this variable. 
       Integer :: nrows_loc, ncols_loc
       Logical :: ludecomposed
       Integer, Pointer :: ipiv (:)
       Complex (8), Pointer :: za (:, :)
#ifdef MPI
       Integer, Dimension(9) :: desc
#endif

    End Type HermitianMatrix
    !
    Type evsystem
       Type (HermitianMatrix) :: hamilton, overlap
    End Type evsystem

    Interface newmatrix
      Module Procedure newComplexMatrix
      Module Procedure newHermitianMatrix
    End Interface

    Interface deletematrix
      Module Procedure deleteComplexMatrix
      Module Procedure deleteHermitianMatrix
    End Interface

    ! temporary interface for the transition to the new function
    Interface HermitianMatrixMatrix
      Module Procedure HermitianMatrixMatrixNew
      Module Procedure HermitianMatrixMatrixOld
    End Interface
    !
  Contains
    !
    !
#ifdef MPI
    Subroutine newComplexMatrix (self, nrows, ncols, distribute)
#else
    Subroutine newComplexMatrix (self, nrows, ncols)
#endif
      Implicit None
      type (ComplexMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: nrows, ncols
#ifdef MPI
      Integer, Intent (In) :: distribute
#endif

      Integer, External   :: NUMROC

!      Integer :: distribute

      self%nrows = nrows
      self%ncols = ncols
#ifdef MPI
!       if (present(distr)) then
!        distribute = distr
! write (*,*) 'distr specified!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!       else
!         distribute = DISTRIBUTE_BOTH
!       end if



      if (distribute .eq. DISTRIBUTE_2D) then
        self%nrows_loc = NUMROC(self%nrows, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
        self%ncols_loc = NUMROC(self%ncols, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
        CALL DESCINIT(self%desc, self%nrows, self%ncols, &
                      MPIglobal%blocksize, MPIglobal%blocksize, 0, 0, &
                      MPIglobal%context, self%nrows_loc, MPIglobal%ierr)
 
      else if (distribute .eq. DISTRIBUTE_ROWS) then
        self%nrows_loc = NUMROC(self%nrows, MPIglobal_1D%blocksize, MPIglobal_1D%myprocrow, 0, MPIglobal_1D%nprocrows)
        self%ncols_loc = self%ncols
        CALL DESCINIT(self%desc, self%nrows, self%ncols, &
                      MPIglobal_1D%blocksize, MPIglobal_1D%blocksize, 0, 0, &
                      MPIglobal_1D%context, self%nrows_loc, MPIglobal_1D%ierr)

       else if (distribute .eq. DISTRIBUTE_COLS) then
         self%nrows_loc = self%nrows
         self%ncols_loc = NUMROC(self%ncols, MPIglobal_1D%blocksize, MPIglobal_1D%myproccol, 0, MPIglobal_1D%nproccols)
        CALL DESCINIT(self%desc, self%nrows, self%ncols, &
                      self%nrows, MPIglobal_1D%blocksize, 0, 0, &
                      MPIglobal_1D%context, self%nrows_loc, MPIglobal_1D%ierr)
      end if 

!       CALL DESCINIT(self%desc, self%nrows, self%ncols, &
!                     MPIglobal%blocksize, MPIglobal%blocksize, 0, 0, &
!                     MPIglobal%context, self%nrows_loc, MPIglobal%ierr)
#else
      self%nrows_loc = nrows
      self%ncols_loc = ncols
#endif
      Allocate (self%za(self%nrows_loc, self%ncols_loc))
      self%za = 0.0

    End Subroutine newComplexMatrix
    !
    !
    Subroutine newHermitianMatrix (self,  size)
!    Subroutine newMatrix (self,  size)
      Implicit None
      type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: size

      Integer, External   :: NUMROC

      self%size = size
#ifdef MPI
      self%nrows_loc = NUMROC(self%size, MPIglobal%blocksize, MPIglobal%myprocrow, 0, MPIglobal%nprocrows)
      self%ncols_loc = NUMROC(self%size, MPIglobal%blocksize, MPIglobal%myproccol, 0, MPIglobal%nproccols)
      CALL DESCINIT(self%desc, self%size, self%size, &
                    MPIglobal%blocksize, MPIglobal%blocksize, 0, 0, &
                    MPIglobal%context, self%nrows_loc, MPIglobal%ierr)
#else
      self%nrows_loc = size
      self%ncols_loc = size
#endif
      self%ludecomposed = .False.

      Allocate (self%za(self%nrows_loc, self%ncols_loc))
      self%za = 0.0

!    End Subroutine newMatrix
    End Subroutine newHermitianMatrix
    !
    !
    Subroutine deleteHermitianMatrix (self)
      Type (HermitianMatrix), Intent (Inout) :: self

      Deallocate (self%za)

      If (self%ludecomposed) deallocate (self%ipiv)
    End Subroutine deleteHermitianMatrix
    !
    !
    Subroutine deleteComplexMatrix (self)
      Type (ComplexMatrix), Intent (Inout) :: self

      Deallocate (self%za)

    End Subroutine deleteComplexMatrix
    !
    !
    Subroutine newsystem (self,   size)
      Type (evsystem), Intent (Out) :: self

      Integer, Intent (In) :: size
      Call newmatrix (self%hamilton, size)
      Call newmatrix (self%overlap,  size)
    End Subroutine newsystem
    !
    !
    Subroutine deletesystem (self)
      Type (evsystem), Intent (Inout) :: self
      Call deletematrix (self%hamilton)
      Call deletematrix (self%overlap)
    End Subroutine deletesystem
    !
    !
! never used???
!     Subroutine Hermitianmatrix_size2update (self, n, alpha, x, y)
!       Type (HermitianMatrix), Intent (Inout) :: self
!       Integer, Intent (In) :: n
!       Complex (8), Intent (In) :: alpha, x (:), y (:)
!       !
! 
!       Call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%size)
! 
!     End Subroutine Hermitianmatrix_size2update
    !
    !
    Subroutine Hermitianmatrix_indexedupdate (self, i, j, z)
      Type (HermitianMatrix), Intent (Inout) :: self
      Integer :: i, j
      Complex (8) :: z
      Integer :: ipx

      If (j .Le. i) Then
         self%za (j, i) = self%za(j, i) + z
      Else
         Write (*,*) "warning lower part of hamilton updated"
      End If

      Return
    End Subroutine Hermitianmatrix_indexedupdate
    !
    !
    Subroutine Hermitianmatrixvector (self, alpha, vin, beta, vout)
      Implicit None
      Type (HermitianMatrix), Intent (Inout) :: self
      Complex (8), Intent (In) :: alpha, beta
      Complex (8), Intent (Inout) :: vin (:)
      Complex (8), Intent (Inout) :: vout (:)
      !

      Call zhemv ("U", self%size, alpha, self%za, self%size, vin, &
           1, beta, vout, 1)

    End Subroutine Hermitianmatrixvector
    !
    !
    Function getsize (self)
      Integer :: getsize
      Type (HermitianMatrix) :: self
      getsize = self%size
    End Function getsize
    !
    !
    Subroutine HermitianmatrixLU (self)
      Type (HermitianMatrix) :: self
      Integer :: info
      If ( .Not. self%ludecomposed) allocate (self%ipiv(self%size))
      !
      If ( .Not. self%ludecomposed) Then

         Call ZGETRF (self%size, self%size, self%za, self%size, &
              self%ipiv, info)

         If (info .Ne. 0) Then
            Write (*,*) "error in iterativearpacksecequn  Hermitianm&
                 &atrixLU "                , info
            Stop
         End If
         self%ludecomposed = .True.
      End If
    End Subroutine HermitianmatrixLU




    subroutine HermitianMatrixMatrixNew(self,zm1,zm2,ldzm,naa,ngp)
! !INPUT/OUTPUT PARAMETERS:
!   self   : Result matrix, hermitian, dimension(ngp,ngp). 
!            The result of the matrix-matrix multiplication 
!            of the other matrices is ADDED to self
!   zm1    : First factor, complex matrix, dimension(naa,ngp)
!   zm2    : Second factor, complex matrix, dimension(naa,ngp)
!   ldzm   : leading dimension of zm1 and zm2
!   naa    : number of rows of zm1 and zm2
!   ngp    : dimension of self
! !DESCRIPTION:
!   Performs the matrix-matrix multiplication
!   self = zm1**H * zm2 + za
!   matrices have dimensions
!     zm1(naa,ngp)
!     zm2(naa,ngp)
!     self%za(ngp,ngp)
!
! !REVISION HISTORY:
!   Created February 2013 (G. Huhs - BSC)
!EOP
!BOC
      Implicit None
!
! arguments
      Type (HermitianMatrix),intent(inout) :: self
      Type (ComplexMatrix),intent(in) :: zm1,zm2
      integer,intent(in)::ldzm,ngp,naa
!
! local variables
      complex(8)::zone=(1.0,0.0)

#ifdef MPI
        CALL PZGEMM( 'C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                     ngp, &           ! M ... rows of op( A ) = rows of C
                     ngp, &           ! N ... cols of op( B ) = cols of C
                     naa, &           ! K ... cols of op( A ) = rows of op( B )
                     zone, &          ! alpha
                     zm1%za, &        ! A, local part
                     1, &             ! IA,
                     1, &             ! JA,
                     zm1%desc, &      ! DESCA,
                     zm2%za, &        ! B, local part
                     1, &             ! IB,
                     1, &             ! JB,
                     zm2%desc, &      ! DESCB,
                     zone, &          ! beta
                     self%za, &       ! C,
                     1, &             ! IC,
                     1, &             ! JC,
                     self%desc &      ! DESCC
                     )
#else
      ! ZGEMM  performs one of the matrix-matrix operations
      !        C := alpha*op( A )*op( B ) + beta*C,
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 ngp, &           ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 naa, &           ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 zm1%za, &        ! A
                 ldzm,&           ! LDA ... leading dimension of A
                 zm2%za, &        ! B
                 ldzm, &          ! LDB ... leading dimension of B
                 zone, &          ! beta
                 self%za, &       ! C
                 self%size &      ! LDC ... leading dimension of C
                )

#endif

    end subroutine HermitianMatrixMatrixNew

    ! Performs the matrix-matrix multiplication
    ! self%za = zm1**H * zm2 + self%za
    ! matrices have dimensions
    !   zm1(naa,ngp)
    !   zm2(naa,ngp)
    !   self%za(ngp,ngp)
    subroutine HermitianMatrixMatrixOld(self,zm1,zm2,ldzm,naa,ngp)
      Implicit None
      Type (HermitianMatrix),intent(inout) :: self
      Complex(8),intent(in)::zm1(:,:),zm2(:,:)
      integer,intent(in)::ldzm,ngp,naa
     
      complex(8)::zone=(1.0,0.0)

      ! ZGEMM  performs one of the matrix-matrix operations
      !        C := alpha*op( A )*op( B ) + beta*C,
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                 ngp, &           ! M ... rows of op( A ) = rows of C
                 ngp, &           ! N ... cols of op( B ) = cols of C
                 naa, &           ! K ... cols of op( A ) = rows of op( B )
                 zone, &          ! alpha
                 zm1, &           ! A
                 ldzm,&           ! LDA ... leading dimension of A
                 zm2, &           ! B
                 ldzm, &          ! LDB ... leading dimension of B
                 zone, &          ! beta
                 self%za(1,1), &  ! C
                 self%size &      ! LDC ... leading dimension of C
                )

    end subroutine HermitianMatrixMatrixOld

    !
    !
    Subroutine Hermitianmatrixlinsolve (self, b)
      Type (HermitianMatrix) :: self
      Complex (8), Intent (Inout) :: b (:)
      Integer :: info
      If (self%ludecomposed) Then

         Call ZGETRS ('N', self%size, 1, self%za, self%size, &
              self%ipiv, b, self%size, info)

         If (info .Ne. 0) Then
            Write (*,*) "error in iterativearpacksecequn Hermitianma&
                 &trixlinsolve "                , info
            Stop
         End If
      End If
    End Subroutine Hermitianmatrixlinsolve
    !
    !
    Subroutine HermitianMatrixAXPY (alpha, x, y)
      Complex (8) :: alpha
      Type (HermitianMatrix) :: x, y
      Integer :: mysize

      mysize = x%size * (x%size)
      Call zaxpy (mysize, alpha, x%za, 1, y%za, 1)

    End Subroutine HermitianMatrixAXPY
    !
    !
    Subroutine HermitianMatrixcopy (x, y)
      Complex (8) :: alpha
      Type (HermitianMatrix) :: x, y
      Integer :: mysize

      mysize = x%size * (x%size)
      Call zcopy (mysize, x%za, 1, y%za, 1)

    End Subroutine HermitianMatrixcopy
    !
    !
    Subroutine HermitianMatrixToFiles (self, prefix)
      Implicit None
      Type (HermitianMatrix), Intent (In) :: self
      Character (256), Intent (In) :: prefix
      Character (256) :: filename

      filename = trim (prefix) // ".real.OUT"
      Open (888, File=filename)
      Write (888,*) dble (self%za)

      Close (888)
      !

      filename = trim (prefix) // ".imag.OUT"
      Open (888, File=filename)
      Write (888,*) aimag (self%za)

      Close (888)
    End Subroutine HermitianMatrixToFiles
    !
    !
    Subroutine HermitianMatrixTruncate (self, threshold)
      Implicit None
      Type (HermitianMatrix), Intent (Inout) :: self
      Real (8), Intent (In) :: threshold
      Integer :: n, i, j
      n = self%size

      Do j = 1, n
         Do i = 1, n
            If (Abs(dble(self%za(i, j))) .Lt. threshold) &
                 self%za(i, j) = self%za(i, j) - dcmplx &
                 (dble(self%za(i, j)), 0)
            If (Abs(aimag(self%za(i, j))) .Lt. threshold) &
                 self%za(i, j) = self%za(i, j) - dcmplx (0, &
                 aimag(self%za(i, j)))
         End Do
      End Do

    End Subroutine HermitianMatrixTruncate
    !
    !
    Subroutine HermitianMatrixdiagonal (self, d)
      Implicit None
      Type (HermitianMatrix), Intent (In) :: self
      Complex (8), Intent (Out) :: d (self%size)
      Integer :: i

      Do i = 1, self%size
         d (i) = self%za(i, i)
      End Do

    End Subroutine HermitianMatrixdiagonal
    !

  End Module modfvsystem
