!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olpaan (overlap, is, ia, ngp, apwalm)
      Use modmain
      Use modinput
      Use modfvsystem
      Implicit None
! arguments
      Type (hermitianmatrix), Intent (Inout) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
!      Complex (8) :: x (ngp), y (ngp)
!
! local variables
      Integer :: ias, l, m, lm, io,naa
      Complex(8),allocatable::zm1(:,:)
! external functions
!      Complex (8) zdotu
!      External zdotu
      ias = idxas (ia, is)
      naa=0
     allocate(zm1(apwordmax*lmmaxapw,ngp))
	 zm1=zzero
	 naa=0
      Do l = 0, input%groundstate%lmaxmat
         Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
              naa=naa+1
              zm1(naa,:)=apwalm(1:ngp, io, lm, ias)
            End Do
         End Do
      End Do
      call HermitianMatrixMatrix(overlap,zm1,zm1,apwordmax*lmmaxapw,naa,ngp)

     deallocate(zm1)

	 Return
End Subroutine
