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
#ifdef MPI
Subroutine hmlistln (hamilton, ngp, igpig, vgpc, ngp_loc_row, ngp_loc_col, rows_loc2glob, cols_loc2glob)
#else
Subroutine hmlistln (hamilton, ngp, igpig, vgpc)
#endif
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
      Integer,                Intent (In)    :: ngp
      Integer,                Intent (In)    :: igpig (ngkmax)
      Real (8),               Intent (In)    :: vgpc (3, ngkmax)
#ifdef MPI      
      Integer, Intent (In)     :: ngp_loc_row, ngp_loc_col !TODO better name
      Integer, Dimension(hamilton%nrows_loc), Intent (In) :: rows_loc2glob
      Integer, Dimension(hamilton%ncols_loc), Intent (In) :: cols_loc2glob
#endif
!
!
      Complex (8) :: zt
! local variables
      Integer :: i_loc, i_glob, j_loc, j_glob, ig, iv (3)
      Integer :: i, j
      Real (8) :: t1
!
! calculate the matrix elements
!#$omp parallel default(shared) &
!#$omp shared(h) private(iv,ig,t1,i,j)
!#$omp do
#ifdef MPI
      Do j_loc = 1, ngp_loc_col
         j_glob = cols_loc2glob(j_loc)
         Do i_loc = 1, ngp_loc_row
            i_glob = rows_loc2glob(i_loc)
            if (i_glob .Le. j_glob) then 
#else
      Do j_glob = 1, ngp
         Do i_glob = 1, j_glob
#endif
               iv (:) = ivg (:, igpig(i_glob)) - ivg (:, igpig(j_glob))
               ig = ivgig (iv(1), iv(2), iv(3))
               If ((ig .Gt. 0) .And. (ig .Le. ngvec)) Then
                  t1 = 0.5d0 * dot_product (vgpc(:, i_glob), vgpc(:, j_glob))
                  zt = veffig (ig) + t1 * cfunig (ig)
#ifdef MPI
                  hamilton%za(i_loc,j_loc) = hamilton%za(i_loc,j_loc) + zt
#else
                  Call Hermitianmatrix_indexedupdate (hamilton, j_glob, i_glob, zt)
#endif
               End If
#ifdef MPI
            End If
#endif
         End Do
      End Do
! #else
!       Do j = 1, ngp
!          Do i = 1, j
!             iv (:) = ivg (:, igpig(i)) - ivg (:, igpig(j))
!             ig = ivgig (iv(1), iv(2), iv(3))
!             If ((ig .Gt. 0) .And. (ig .Le. ngvec)) Then
!                t1 = 0.5d0 * dot_product (vgpc(:, i), vgpc(:, j))
!                zt = veffig (ig) + t1 * cfunig (ig)
! 
!                Call Hermitianmatrix_indexedupdate (hamilton, j, i, zt)
!             End If
!          End Do
!       End Do
! #endif
!#$omp end do
!#$omp end parallel
!
      Return
End Subroutine
!EOC
