!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlistl
! !INTERFACE:
!
!
Subroutine hmlistln (hamilton, ngp, igpig, vgpc)
! !USES:
      Use modmain
      Use modfvsystem
! !INPUT/OUTPUT PARAMETERS:
!   ngp   : number of G+p-vectors (in,integer)
!   igpig : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc  : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   hamilton : the Hamiltonian matrix in packed form (inout,complex(npmatmax))
! !DESCRIPTION:
!   Computes the interstitial contribution to the Hamiltonian matrix for the APW
!   basis functions. The Hamiltonian is given by
!   $$ H^{\rm I}({\bf G+k,G'+k})=\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    \tilde{\Theta}({\bf G-G'})+V^{\sigma}({\bf G-G'}), $$
!   where $V^{\sigma}$ is the effective interstitial potential for spin $\sigma$
!   and $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Type (HermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
!
!
      Complex (8) :: zt
! local variables
      Integer :: i, j, k, ig, iv (3)
      Real (8) :: t1
      Complex (8) zt1
!
! calculate the matrix elements
!#$omp parallel default(shared) &
!#$omp shared(h) private(iv,ig,t1,i,j)
!#$omp do
      Do j = 1, ngp
        Do i = 1, j
          iv (:) = ivg (:, igpig(i)) - ivg (:, igpig(j))
          ig = ivgig (iv(1), iv(2), iv(3))
          If ((ig .Gt. 0) .And. (ig .Le. ngvec)) Then
            t1 = 0.5d0 * dot_product (vgpc(:, i), vgpc(:, j))
!            write(*,*) t1,vgpc(:, 1)
!            stop
            zt = veffig (ig) + t1 * cfunig (ig)

            Call Hermitianmatrix_indexedupdate (hamilton, j, i, zt)
          End If
        End Do
      End Do
!#$omp end do
!#$omp end parallel
!
      Return
End Subroutine
!EOC
