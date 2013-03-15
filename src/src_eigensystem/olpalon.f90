!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine olpalon (system, is, ia, apwalm)
      Use mod_muffin_tin,  Only: idxlm, lmmaxapw
      Use mod_APW_LO,      Only: apword, apwordmax, lorbl, nlorb
      Use mod_atoms,       Only: natmtot, idxas
      Use mod_eigensystem, Only: oalo, idxlo
      Use modfvsystem,     Only: evsystem,Hermitianmatrix_indexedupdate
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Complex (8), Intent (In) :: apwalm(system%ngp_loc_rows, apwordmax, lmmaxapw, natmtot) !SPLIT first dimension over procs
!
! local variables
      Integer :: ias, ilo, io, l, m, lm, i, j_glob, j_loc
      Complex (8) zsum

      j_loc = system%ngp_loc_cols
      ias = idxas(ia, is)

      Do ilo = 1, nlorb (is)
         l = lorbl (ilo, is)
         Do m = - l, l
            lm = idxlm (l, m)
            j_glob = system%ngp + idxlo (lm, ilo, ias)
#ifdef MPI
            if (Any(j_glob .eq. system%overlap%my_cols_idx)) then
               j_loc = j_loc+1 !!! assuming that the indexing increases
                                 ! => problem if indexing or order of loops changes
                                 ! an array linking between global and local indices would be better
#endif

! calculate the matrix elements
               Do i = 1, system%ngp_loc_rows
                  zsum = 0.d0
                  Do io = 1, apword (l, is)
                     zsum = zsum + Conjg(apwalm(i, io, lm, ias)) * oalo &
                     & (io, ilo, ias)
                  End Do
#ifdef MPI
                  system%overlap%za(i,j_loc) = system%overlap%za(i,j_loc) + zsum
#else
                  Call Hermitianmatrix_indexedupdate (system%overlap, j_glob, i, zsum)
#endif
               End Do
#ifdef MPI
            End If
#endif
         End Do
      End Do

      Return
End Subroutine
