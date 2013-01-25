 
  ! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
  ! This file is distributed under the terms of the GNU General Public License.
  ! See the file COPYING for license details.

  ! Module for setting up the eigensystem
  ! it is designed in a way that all other subroutines
  ! dealing with setting up and solving the system can acsess the
  ! data transparently allowing to choose from different datatypes
  ! more easily
  Module modfvsystem
    Implicit None
    !
    Type HermitianMatrix
       !
       Integer :: size   !TODO: in the end we should be able to get rid of this variable. 
       Integer :: nrows_loc, ncols_loc
       Logical :: ludecomposed
       Integer, Pointer :: ipiv (:)
       Complex (8), Pointer :: za (:, :)

    End Type HermitianMatrix
    !
    Type evsystem
       Type (HermitianMatrix) :: hamilton, overlap
    End Type evsystem
    !
  Contains
    !
    !
    Subroutine newmatrix (self,  size)
#ifdef MPI
      use modmpi
#endif
      type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: size

      Integer, External   :: NUMROC

      self%size = size
#ifdef MPI
        self%nrows_loc = NUMROC(self%size, blocksize, myprocrow, 0, nprocrows)
        self%ncols_loc = NUMROC(self%size, blocksize, myproccol, 0, nproccols)
#else
        self%nrows_loc = size
        self%ncols_loc = size
#endif
      self%ludecomposed = .False.

      Allocate (self%za(self%nrows_loc, self%ncols_loc))
      self%za = 0.0

    End Subroutine newmatrix
    !
    !
    Subroutine deletematrix (self)
      Type (HermitianMatrix), Intent (Inout) :: self

      Deallocate (self%za)

      If (self%ludecomposed) deallocate (self%ipiv)
    End Subroutine deletematrix
    !
    !
    Subroutine newsystem (self,   size)
      Type (evsystem), Intent (Out) :: self

      Integer, Intent (In) :: size
      Call newmatrix (self%hamilton,   size)
      Call newmatrix (self%overlap,  size)
    End Subroutine newsystem
    !
    !
    Subroutine deleteystem (self)
      Type (evsystem), Intent (Inout) :: self
      Call deletematrix (self%hamilton)
      Call deletematrix (self%overlap)
    End Subroutine deleteystem
    !
    !
    Subroutine Hermitianmatrix_size2update (self, n, alpha, x, y)
      Type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: alpha, x (:), y (:)
      !

      Call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%size)

    End Subroutine Hermitianmatrix_size2update
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

    subroutine HermitianMatrixMatrix(self,zm1,zm2,ldzm,naa,ngp)
     Type (HermitianMatrix),intent(inout) :: self
     Complex(8),intent(in)::zm1(:,:),zm2(:,:)
     integer,intent(in)::ldzm,ngp,naa
     complex(8)::zone=(1.0,0.0)
     call zgemm('C','N',ngp,ngp,naa,zone,zm1,ldzm,&
  	   			 zm2,ldzm,zone,self%za(1,1),self%size)
    end subroutine

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
