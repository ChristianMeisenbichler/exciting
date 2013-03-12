!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
#ifdef MPI      
Subroutine olpalon (overlap, is, ia, ngp, apwalm, ngp_loc_row, ngp_loc_col, overlap_loc_cols)
#else
Subroutine olpalon (overlap, is, ia, ngp, apwalm)
#endif
      Use modmain
      Use modfvsystem
      Implicit None
! arguments
      Type (HermitianMatrix) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
#ifdef MPI      
      Complex (8), Intent (In) :: apwalm (ngp_loc_row, apwordmax, lmmaxapw, natmtot) !SPLIT first dimension over procs
      Integer, Intent (In)     :: ngp_loc_row, ngp_loc_col !TODO better name
      Integer, Dimension(overlap%ncols_loc), Intent (In) :: overlap_loc_cols
#else 
      Complex (8), Intent (In) :: apwalm (ngp, apwordmax, lmmaxapw, natmtot) 
#endif
!
! local variables
      Integer :: ias, ilo, io, l, m, lm, i, j_glob, j_loc, max_i
      Complex (8) zsum

#ifdef MPI
      max_i = ngp_loc_row
      j_loc = ngp_loc_col
#else
      max_i = ngp
#endif

      ias = idxas (ia, is)
      Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         Do m = - l, l
            lm = idxlm (l, m)
            j_glob = ngp + idxlo (lm, ilo, ias)
#ifdef MPI
            if (Any(overlap_loc_cols .eq. j_glob)) then
               j_loc = j_loc+1 !!! assuming that the indexing increases
                                 ! => problem if indexing or order of loops changes
                                 ! an array linking between global and local indices would be better
#endif

! calculate the matrix elements
               Do i = 1, max_i
                  zsum = 0.d0
                  Do io = 1, apword (l, is)
                     zsum = zsum + conjg (apwalm(i, io, lm, ias)) * oalo &
                     & (io, ilo, ias)
                  End Do
#ifdef MPI
                  overlap%za(i,j_loc) = overlap%za(i,j_loc) + zsum
#else
                  Call Hermitianmatrix_indexedupdate (overlap, j_glob, i, zsum)
#endif
               End Do
#ifdef MPI
            End If
#endif
         End Do
      End Do

      Return
End Subroutine
