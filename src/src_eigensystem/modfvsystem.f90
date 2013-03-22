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
!   Parallelized, February 2013 (G. Huhs - BSC)
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
       Complex (8), Pointer :: za (:,:)
#ifdef MPI
       Integer, Dimension(:), Pointer :: desc
       Integer :: distribute
#endif
    End Type ComplexMatrix

    Type HermitianMatrix
       !
       Integer :: size   !TODO: in the end we should be able to get rid of this variable. 
       Integer :: nrows_loc, ncols_loc
       Logical :: ludecomposed
       Integer, Pointer :: ipiv (:)
       Complex (8), Pointer :: za (:,:)
#ifdef MPI
       Integer, Dimension(:), Pointer :: desc
       Integer, Dimension(:), Pointer :: my_rows_idx, my_cols_idx
       Integer :: distribute
#endif
    End Type HermitianMatrix
    !
    Type evsystem
       Type (HermitianMatrix) :: hamilton, overlap

       ! number of G+k-vectors for augmented plane waves
       Integer                :: ngp
       ! number of local rows/cols of G+k vectors
       Integer                :: ngp_loc_rows, ngp_loc_cols

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
      Integer, Intent (In), Optional :: distribute

      Integer, External   :: NUMROC
#endif

      self%nrows = nrows
      self%ncols = ncols
#ifdef MPI
      if (present(distribute)) then
       self%distribute = distribute
      else
        self%distribute = DISTRIBUTE_2D
      end if

      allocate(self%desc(9))

      select case (self%distribute)
        case (DISTRIBUTE_2D) 
          self%nrows_loc = NUMROC(self%nrows, MPIkpt_2D%blocksize, MPIkpt_2D%myprocrow, 0, MPIkpt_2D%nprocrows)
          self%ncols_loc = NUMROC(self%ncols, MPIkpt_2D%blocksize, MPIkpt_2D%myproccol, 0, MPIkpt_2D%nproccols)
          CALL DESCINIT(self%desc, self%nrows, self%ncols, &
                        MPIkpt_2D%blocksize, MPIkpt_2D%blocksize, 0, 0, &
                        MPIkpt_2D%context, self%nrows_loc, MPIkpt_2D%ierr)
  
!         case (DISTRIBUTE_ROWS)
!           self%nrows_loc = NUMROC(self%nrows, MPIkpt_1D%blocksize, MPIkpt_1D%myprocrow, 0, MPIkpt_1D%nprocrows)
!           self%ncols_loc = self%ncols
!           CALL DESCINIT(self%desc, self%nrows, self%ncols, &
!                         MPIkpt_1D%blocksize, self%ncols, 0, 0, &
!                         MPIkpt_1D%context, self%nrows_loc, MPIkpt_1D%ierr)

        case (DISTRIBUTE_COLS) 
          self%nrows_loc = self%nrows
          self%ncols_loc = NUMROC(self%ncols, MPIkpt_1D%blocksize, MPIkpt_1D%myproccol, 0, MPIkpt_1D%nproccols)
          CALL DESCINIT(self%desc, self%nrows, self%ncols, &
                        self%nrows, MPIkpt_1D%blocksize, 0, 0, &
                        MPIkpt_1D%context, self%nrows_loc, MPIkpt_1D%ierr)
      end select 

#else
      self%nrows_loc = nrows
      self%ncols_loc = ncols
#endif
      Allocate (self%za(self%nrows_loc, self%ncols_loc))
      self%za = cmplx(0,0,8)

    End Subroutine newComplexMatrix
    !
    !
#ifdef MPI
    Subroutine newHermitianMatrix (self, size, distribute)
#else
    Subroutine newHermitianMatrix (self,  size)
#endif
      Implicit None
      type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: size
#ifdef MPI
      Integer, Intent (In), Optional :: distribute

      Integer :: i

      Integer, External   :: NUMROC
#endif

      self%size = size
#ifdef MPI
      if (present(distribute)) then
        self%distribute = distribute
      else
        self%distribute = DISTRIBUTE_2D
      end if

      allocate(self%desc(9))

      select case (self%distribute)
        case (DISTRIBUTE_2D) 
          self%nrows_loc = NUMROC(self%size, MPIkpt_2D%blocksize, MPIkpt_2D%myprocrow, 0, MPIkpt_2D%nprocrows)
          self%ncols_loc = NUMROC(self%size, MPIkpt_2D%blocksize, MPIkpt_2D%myproccol, 0, MPIkpt_2D%nproccols)
          CALL DESCINIT(self%desc, self%size, self%size, &
                        MPIkpt_2D%blocksize, MPIkpt_2D%blocksize, 0, 0, &
                        MPIkpt_2D%context, self%nrows_loc, MPIkpt_2D%ierr)

          Allocate(self%my_rows_idx(self%nrows_loc), self%my_cols_idx(self%ncols_loc))
          Call getLocalIndices(self%my_rows_idx, self%size, MPIkpt_2D%blocksize, MPIkpt_2D%myprocrow, MPIkpt_2D%nprocrows)
          Call getLocalIndices(self%my_cols_idx, self%size, MPIkpt_2D%blocksize, MPIkpt_2D%myproccol, MPIkpt_2D%nproccols)

!         case (DISTRIBUTE_ROWS)
!           self%nrows_loc = NUMROC(self%size, MPIkpt_1D%blocksize, MPIkpt_1D%myprocrow, 0, MPIkpt_1D%nprocrows)
!           self%ncols_loc = self%size
!           CALL DESCINIT(self%desc, self%size, self%size, &
!                         MPIkpt_1D%blocksize, self%size, 0, 0, &
!                         MPIkpt_1D%context, self%nrows_loc, MPIkpt_1D%ierr)
! 
!           Allocate(self%my_rows_idx(self%nrows_loc), self%my_cols_idx(self%ncols_loc))
!           Call getLocalIndices(self%my_rows_idx, self%size, MPIkpt_1D%blocksize, MPIkpt_1D%myprocrow, MPIkpt_1D%nprocrows)
!           self%my_cols_idx = (/(i,i=1,self%ncols_loc)/)

        case (DISTRIBUTE_COLS) 
          self%nrows_loc = self%size
          self%ncols_loc = NUMROC(self%size, MPIkpt_1D%blocksize, MPIkpt_1D%myproccol, 0, MPIkpt_1D%nproccols)
          CALL DESCINIT(self%desc, self%size, self%size, &
                        self%size, MPIkpt_1D%blocksize, 0, 0, &
                        MPIkpt_1D%context, self%nrows_loc, MPIkpt_1D%ierr)

          Allocate(self%my_rows_idx(self%nrows_loc), self%my_cols_idx(self%ncols_loc))
          self%my_rows_idx = (/(i,i=1,self%nrows_loc)/)
          Call getLocalIndices(self%my_cols_idx, self%size, MPIkpt_1D%blocksize, MPIkpt_1D%myproccol, MPIkpt_1D%nproccols)

      end select 
      
#else
      self%nrows_loc = size
      self%ncols_loc = size
#endif
      self%ludecomposed = .False.

      Allocate (self%za(self%nrows_loc, self%ncols_loc))
      self%za = cmplx(0,0,8)

!    End Subroutine newMatrix
    End Subroutine newHermitianMatrix
    !
    !
    Subroutine deleteHermitianMatrix (self)
      Type (HermitianMatrix), Intent (Inout) :: self

      Deallocate (self%za)
      If (self%ludecomposed) deallocate (self%ipiv)

#ifdef MPI
      Deallocate (self%desc)
      Deallocate (self%my_rows_idx, self%my_cols_idx)
#endif
    End Subroutine deleteHermitianMatrix
    !
    !
    Subroutine deleteComplexMatrix (self)
      Type (ComplexMatrix), Intent (Inout) :: self

      Deallocate (self%za)
#ifdef MPI
      Deallocate (self%desc)
#endif
    End Subroutine deleteComplexMatrix
    !
    !
#ifdef MPI
    Subroutine newsystem (self, size, ngp, distribute_prm)
#else
    Subroutine newsystem (self, size, ngp)
#endif
      Type (evsystem), Intent (Out) :: self
      Integer, Intent (In)          :: size, ngp
#ifdef MPI
      Integer, Intent (In), Optional :: distribute_prm

      Integer             :: distribute

      Integer, External   :: NUMROC
#endif

      self%ngp = ngp
#ifdef MPI
      if (present(distribute_prm)) then
        distribute = distribute_prm
      else
        distribute = DISTRIBUTE_2D
      end if

      select case (distribute)
        case (DISTRIBUTE_2D) 
          self%ngp_loc_rows = NUMROC(ngp, MPIkpt_2D%blocksize, MPIkpt_2D%myprocrow, 0, MPIkpt_2D%nprocrows)
          self%ngp_loc_cols = NUMROC(ngp, MPIkpt_2D%blocksize, MPIkpt_2D%myproccol, 0, MPIkpt_2D%nproccols)
!         case (DISTRIBUTE_ROWS)
!           self%ngp_loc_rows = NUMROC(ngp, MPIkpt_1D%blocksize, MPIkpt_1D%myprocrow, 0, MPIkpt_1D%nprocrows)
!           self%ngp_loc_cols = ngp
        case (DISTRIBUTE_COLS) 
          self%ngp_loc_rows = ngp
          self%ngp_loc_cols = NUMROC(ngp, MPIkpt_1D%blocksize, MPIkpt_1D%myproccol, 0, MPIkpt_1D%nproccols)
      end select 
      Call newmatrix (self%hamilton, size, distribute)
      Call newmatrix (self%overlap,  size, distribute)
#else
      self%ngp_loc_rows = ngp
      self%ngp_loc_cols = ngp
      Call newmatrix (self%hamilton, size)
      Call newmatrix (self%overlap,  size)
#endif
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
    Subroutine Hermitianmatrix_indexedupdate (self, col, row, z)
      Type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent(In), Target            :: col, row
      Complex(8), Intent(In)                 :: z

#ifdef MPI
      Integer, Dimension(:), Pointer :: rows_loc2glob
      Integer, Dimension(:), Pointer :: cols_loc2glob
      Integer:: row_glob, col_glob
#else
      Integer, Pointer :: row_glob, col_glob
#endif
     
! Performance in MPI-Mode???
! indexing for each element might cost too much, maybe it's better to do without the checking

#ifdef MPI
      rows_loc2glob => self%my_rows_idx
      cols_loc2glob => self%my_cols_idx
      row_glob = rows_loc2glob(row)
      col_glob = cols_loc2glob(col)
#else
      row_glob => row
      col_glob => col
#endif

      If (row_glob .Le. col_glob) Then
         self%za(row, col) = self%za(row, col) + z
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
!   self   : Result matrix, hermitian. 
!            The result of the matrix-matrix multiplication 
!            of the other matrices is ADDED to self
!   zm1    : First factor, complex matrix, dimension(naa,ngp)
!   zm2    : Second factor, complex matrix, dimension(naa,ngp)
!   ldzm   : leading dimension of zm1 and zm2
!   naa    : number of rows of zm1 and zm2
!   ngp    : number of G+p-vectors, smaller than dimension of self!!
! !DESCRIPTION:
!   Performs the matrix-matrix multiplication
!   self = zm1**H * zm2 + self
!   The matrices are expected to have distributed columns when MPI is used
!      (even for only one proc)
!      Practically just featuring the same distribution (context) is sufficient, 
!      but this is neither needed nor tested
!
! !REVISION HISTORY:
!   Parallelized and with new matrix types February 2013 (G. Huhs - BSC)
!   Created ?? (CHM)
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
      complex(8), dimension(:,:), pointer :: A, B
      integer, dimension(:), pointer      :: descA, descB
!       Type (ComplexMatrix)                :: tmp1, tmp2
!       logical                             :: redistributeA, redistributeB

! PBLAS can't deal with differing processor configurations (contexts)
! so for differing distributions a redistribution is necessary
! here is some code to bring all to 2D distribution assuming a 2D distributed Hamiltonian/overlap matrix 
!       redistributeA = (zm1%distribute .ne. DISTRIBUTE_2D)
!       redistributeB = (zm2%distribute .ne. DISTRIBUTE_2D)
! 
!       if (redistributeA) then 
!         Call newmatrix(tmp1, zm1%nrows, zm1%ncols, DISTRIBUTE_2D)
!         Call PZGEMR2D(zm1%nrows, zm1%ncols, zm1%za, 1, 1, zm1%desc, tmp1%za, 1, 1, tmp1%desc, MPIkpt_2D%context)
!         A => tmp1%za
!         descA => tmp1%desc
!       else 
!         A => zm1%za
!         descA => zm1%desc
!       end if
! 
!       if (redistributeB) then 
!         Call newmatrix(tmp2, zm2%nrows, zm2%ncols, DISTRIBUTE_2D)
!         Call PZGEMR2D(zm2%nrows, zm2%ncols, zm2%za, 1, 1, zm2%desc, tmp2%za, 1, 1, tmp2%desc, MPIkpt_2D%context)
!         B => tmp2%za
!         descB => tmp2%desc
!       else 
!         B => zm2%za
!         descB => zm2%desc
!       end if

      A => zm1%za
      descA => zm1%desc
      B => zm2%za
      descB => zm2%desc

      CALL PZGEMM( 'C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                   ngp, &           ! M ... rows of op( A ) = rows of C
                   ngp, &           ! N ... cols of op( B ) = cols of C
                   naa, &           ! K ... cols of op( A ) = rows of op( B )
                   zone, &          ! alpha
                   A, &             ! A, local part
                   1, &             ! IA,
                   1, &             ! JA,
                   descA, &         ! DESCA,
                   B, &             ! B, local part
                   1, &             ! IB,
                   1, &             ! JB,
                   descB, &         ! DESCB,
                   zone, &          ! beta
                   self%za, &       ! C,
                   1, &             ! IC,
                   1, &             ! JC,
                   self%desc &      ! DESCC
                   )

!       if (redistributeA) then 
!         Call deletematrix(tmp1)
!       end if
!       if (redistributeB) then 
!         Call deletematrix(tmp2)
!       end if
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

    Subroutine GetLocalIndices(loc_idx, n_glob, blocksize, my_procidx, n_procs)
      Implicit None
      Integer, Dimension(:), Intent(Out) :: loc_idx
      Integer, Intent(In) :: n_glob, blocksize, my_procidx, n_procs

      Integer :: i_glob, i_loc

      If (n_procs .eq. 1) Then
        loc_idx = (/(i_glob,i_glob=1,n_glob)/)
      Else
        i_loc = 1
        Do i_glob=1,n_glob
          If (Mod((i_glob-1)/blocksize,n_procs) .eq. my_procidx) Then
            loc_idx(i_loc) = i_glob
            i_loc = i_loc+1
          End If
        End Do
      End If
    End Subroutine GetLocalIndices


    Function glob2loc(loc_indices,glob)
      Implicit None
! return value
      Integer :: glob2loc
! arguments
      Integer, Dimension(:), Intent(In)  :: loc_indices
      Integer,               Intent(In)  :: glob
!
! local variables
      Integer :: num_loc_indices

      num_loc_indices = size(loc_indices)

      Do glob2loc=1,num_loc_indices
         If (loc_indices(glob2loc) .eq. glob) Return
      End Do

    End Function glob2loc

    
#ifdef MPI
    Subroutine RedistributeHermitianMatrix1DTo2D(self)
      Implicit None
! arguments
      Type(HermitianMatrix), Intent(InOut) :: self
! local variables
      Type(HermitianMatrix) :: mat2D

      Call newMatrix(mat2D, self%size, DISTRIBUTE_2D)
      Call PZGEMR2D(self%size, self%size, self%za, 1, 1, self%desc, mat2D%za, 1, 1, mat2D%desc, MPIkpt_2D%context)
      
      Deallocate(self%desc, self%za, self%my_rows_idx, self%my_cols_idx)
      self%desc => mat2D%desc
      self%za => mat2D%za
      self%my_rows_idx => mat2D%my_rows_idx
      self%my_cols_idx => mat2D%my_cols_idx

     self%nrows_loc = mat2D%nrows_loc
     self%ncols_loc = mat2D%ncols_loc

    End Subroutine RedistributeHermitianMatrix1DTo2D
#endif


  End Module modfvsystem
