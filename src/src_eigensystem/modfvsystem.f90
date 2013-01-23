 
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
       Integer :: rank
       Logical ::   ludecomposed
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
    Subroutine newmatrix (self,   rank)
      type (HermitianMatrix), Intent (Inout) :: self

      Integer, Intent (In) :: rank
      self%rank = rank

      self%ludecomposed = .False.

      Allocate (self%za(rank, rank))
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
    Subroutine newsystem (self,   rank)
      Type (evsystem), Intent (Out) :: self

      Integer, Intent (In) :: rank
      Call newmatrix (self%hamilton,   rank)
      Call newmatrix (self%overlap,  rank)
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
    Subroutine Hermitianmatrix_rank2update (self, n, alpha, x, y)
      Type (HermitianMatrix), Intent (Inout) :: self
      Integer, Intent (In) :: n
      Complex (8), Intent (In) :: alpha, x (:), y (:)
      !

      Call ZHER2 ('U', n, alpha, x, 1, y, 1, self%za, self%rank)

    End Subroutine Hermitianmatrix_rank2update
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

      Call zhemv ("U", self%rank, alpha, self%za, self%rank, vin, &
           1, beta, vout, 1)

    End Subroutine Hermitianmatrixvector
    !
    !
    Function getrank (self)
      Integer :: getrank
      Type (HermitianMatrix) :: self
      getrank = self%rank
    End Function getrank
    !
    !
    Subroutine HermitianmatrixLU (self)
      Type (HermitianMatrix) :: self
      Integer :: info
      If ( .Not. self%ludecomposed) allocate (self%ipiv(self%rank))
      !
      If ( .Not. self%ludecomposed) Then

         Call ZGETRF (self%rank, self%rank, self%za, self%rank, &
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
  	   			 zm2,ldzm,zone,self%za(1,1),self%rank)
    end subroutine

    !
    !
    Subroutine Hermitianmatrixlinsolve (self, b)
      Type (HermitianMatrix) :: self
      Complex (8), Intent (Inout) :: b (:)
      Integer :: info
      If (self%ludecomposed) Then

         Call ZGETRS ('N', self%rank, 1, self%za, self%rank, &
              self%ipiv, b, self%rank, info)

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

      mysize = x%rank * (x%rank)
      Call zaxpy (mysize, alpha, x%za, 1, y%za, 1)

    End Subroutine HermitianMatrixAXPY
    !
    !
    Subroutine HermitianMatrixcopy (x, y)
      Complex (8) :: alpha
      Type (HermitianMatrix) :: x, y
      Integer :: mysize

      mysize = x%rank * (x%rank)
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
      n = self%rank

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
      Complex (8), Intent (Out) :: d (self%rank)
      Integer :: i

      Do i = 1, self%rank
         d (i) = self%za(i, i)
      End Do

    End Subroutine HermitianMatrixdiagonal
    !
    Subroutine solvewithlapack(system,nstfv,evecfv,evalfv)
      use mod_timing
      use modinput
      use mod_eigensystem, only:nmatmax
      type(evsystem)::system
      integer::nstfv

      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
      !local
      Integer :: is, ia, i, m, np, info ,nmatp
      Real (8) :: vl, vu
      Real (8) :: ts0, ts1
      ! allocatable arrays
      Integer, Allocatable :: iwork (:)
      Integer, Allocatable :: ifail (:)
      Real (8), Allocatable :: w (:)
      Real (8), Allocatable :: rwork (:)
      Complex (8), Allocatable :: v (:)
      Complex (8), Allocatable :: work (:)
      Call timesec (ts0)

         vl = 0.d0
         vu = 0.d0
         ! LAPACK 3.0 call
         !nmatmax
         nmatp=system%hamilton%rank
         Allocate (iwork(5*nmatp))
         Allocate (ifail(nmatp))
         Allocate (w(nmatp))
         Allocate (rwork(7*nmatp))
         Allocate (v(1))
         Allocate (work(2*nmatp))
         !Call zhpgvx (1, 'V', 'I', 'U', nmatp, system%hamilton%zap, &
         !system%overlap%zap, vl, vu, 1, nstfv, &
         !input%groundstate%solver%evaltol, m, w, evecfv, nmatmax, work, &
         !rwork, iwork, ifail, info)
         call ZHEGVX(1, 'V', 'I', 'U', nmatp, system%hamilton%za, nmatp, system%overlap%za, nmatp,&
              vl, vu, 1, nstfv, input%groundstate%solver%evaltol, &
              m, w, evecfv, nmatmax, work, 2*nmatp, rwork, iwork, ifail, info )
         evalfv (1:nstfv) = w (1:nstfv)
         !
         !
         !
         If (info .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(seceqnfv): diagonalisation failed")')
            Write (*, '(" ZHPGVX returned INFO = ", I8)') info
            If (info .Gt. nmatp) Then
               i = info - nmatp
               Write (*, '(" The leading minor of the overlap matrix of or&
                    &der ", I8)'                ) i
               Write (*, '("  is not positive definite")')
               Write (*, '(" Order of overlap matrix : ", I8)') nmatp
               Write (*,*)
            End If
            Stop
         End If
!         Call timesec (ts1)
         !$OMP CRITICAL
!         timefv = timefv + ts1 - ts0
         !$OMP END CRITICAL
         Call deleteystem (system)
         Deallocate (iwork, ifail, w, rwork, v, work)
         Call timesec (ts1)
         timefv = timefv + ts1 - ts0       
    end subroutine solvewithlapack
  End Module modfvsystem
