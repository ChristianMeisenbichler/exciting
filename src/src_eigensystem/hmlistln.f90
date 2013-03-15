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
Subroutine hmlistln (system, igpig, vgpc)
! !USES:
    Use mod_Gkvector,    Only: ngkmax
    Use mod_Gvector,     Only: ivg,ngvec,cfunig,ivgig
    Use mod_potential_and_density, Only : veffig
    Use modfvsystem,     Only: evsystem, Hermitianmatrix_indexedupdate

! !INPUT/OUTPUT PARAMETERS:
!   system : EVSystem with the Hamiltonian to update (inout,evsystem)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!
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
!   Parallelized and with new interface March 2013 (G. Huhs - BSC)
!   Created April 2003 (JKD)

!EOP
!BOC
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer,         Intent (In)    :: igpig (ngkmax)
      Real(8),        Intent (In)     :: vgpc (3, ngkmax)
!
! local variables
#ifdef MPI
      Integer, Dimension(:), Pointer :: rows_loc2glob
      Integer, Dimension(:), Pointer :: cols_loc2glob
      Integer :: i_loc, j_loc
      Integer :: i_glob, j_glob
#else
      Integer, Pointer :: i_loc, j_loc
      Integer, Target  :: i_glob, j_glob
#endif
      Complex(8) :: zt
      Integer    :: ig, iv(3)
      Real(8)    :: t1
!
! calculate the matrix elements
!#$omp parallel default(shared) &
!#$omp shared(h) private(iv,ig,t1,i,j)
!#$omp do
#ifdef MPI
      rows_loc2glob => system%hamilton%my_rows_idx
      cols_loc2glob => system%hamilton%my_cols_idx
      Do j_loc = 1, system%ngp_loc_cols
         j_glob = cols_loc2glob(j_loc)
         Do i_loc = 1, system%ngp_loc_rows
            i_glob = rows_loc2glob(i_loc)
            if (i_glob .Le. j_glob) then 
#else
      i_loc => i_glob
      j_loc => j_glob
      Do j_glob = 1, system%ngp
         Do i_glob = 1, j_glob
#endif
               iv (:) = ivg (:, igpig(i_glob)) - ivg (:, igpig(j_glob))
               ig = ivgig (iv(1), iv(2), iv(3))
               If ((ig .Gt. 0) .And. (ig .Le. ngvec)) Then
                  t1 = 0.5d0 * dot_product (vgpc(:, i_glob), vgpc(:, j_glob))
                  zt = veffig (ig) + t1 * cfunig (ig)
#ifdef MPI
                  system%hamilton%za(i_loc,j_loc) = system%hamilton%za(i_loc,j_loc) + zt
#else
                  Call Hermitianmatrix_indexedupdate(system%hamilton, j_loc, i_loc, zt)
#endif
               End If
#ifdef MPI
            Else
               Exit
            End If
#endif
         End Do
      End Do
!#$omp end do
!#$omp end parallel
!
      Return
End Subroutine
!EOC
