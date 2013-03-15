!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: olpistl
! !INTERFACE:
!
!
Subroutine olpistln (system, igpig)
! !USES:
     Use mod_Gkvector,    Only: ngkmax
     Use mod_Gvector,     Only: ivg,ngvec,cfunig,ivgig
     Use modfvsystem,     Only: evsystem,Hermitianmatrix_indexedupdate
! !INPUT/OUTPUT PARAMETERS:
!   system : EVsystem with the overlap matrix to update (inout,evsystem)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
! !DESCRIPTION:
!   Computes the interstitial contribution to the overlap matrix for the APW
!   basis functions. The overlap is given by
!   $$ O^{\rm I}({\bf G+k,G'+k})=\tilde{\Theta}({\bf G-G'}), $$
!   where $\tilde{\Theta}$ is the characteristic function. See routine
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Parallelized and new interface March 2013 (Georg Huhs, BSC)
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer, Intent (In) :: igpig (ngkmax)
!
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
      Integer :: iv(3), ig
!
! calculate the matrix elements
!$omp parallel default(shared) &
!$omp  private(iv,ig,i,j)
!$omp do
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
               Call Hermitianmatrix_indexedupdate(system%overlap, j_loc, i_loc, cfunig(ig))
!commented: performance-version for MPI
! #ifdef MPI
!                system%overlap%za(i_loc,j_loc) = system%overlap%za(i_loc,j_loc) + cfunig(ig)
! #else
!                Call Hermitianmatrix_indexedupdate(system%overlap, j_glob, i_glob, cfunig(ig))
! #endif
               End If
#ifdef MPI
            Else
               Exit
            End If
#endif
        End Do
      End Do
!$omp end do
!$omp end parallel
!
      Return
End Subroutine
!EOC
