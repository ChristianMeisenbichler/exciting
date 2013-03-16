!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!

Subroutine olplolon (system, is, ia)
#ifdef MPI
use modmpi
#endif
      Use mod_muffin_tin,  Only: idxlm
      Use mod_APW_LO,      Only: lorbl, nlorb
      Use mod_atoms,       Only: idxas
      Use mod_eigensystem, Only: ololo, idxlo
      Use modfvsystem,     Only: evsystem, Hermitianmatrix_indexedupdate,glob2loc
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
!
!
! local variables
#ifdef MPI
      Integer :: i_loc, j_loc
      Integer :: i_glob, j_glob
#else
      Integer, Pointer :: i_loc, j_loc
      Integer, Target  :: i_glob, j_glob
#endif
      Complex (8) :: zt
      Integer :: ias, ilo1, ilo2, l, m, lm
      ias = idxas (ia, is)

#ifndef MPI
      i_loc => i_glob
      j_loc => j_glob
#endif
#ifdef MPI
write (*,*) MPIglobal%rank, 'loc2glob', system%overlap%my_rows_idx
#endif
      Do ilo1 = 1, nlorb (is)
         l = lorbl (ilo1, is)
         Do ilo2 = 1, nlorb (is)
            If (lorbl(ilo2, is) .Eq. l) Then
               Do m = - l, l
                  lm = idxlm (l, m)
                  i_glob = system%ngp + idxlo (lm, ilo1, ias)
                  j_glob = system%ngp + idxlo (lm, ilo2, ias)
                  If (i_glob .Le. j_glob) Then
#ifdef MPI
                     If (Any(i_glob .eq. system%overlap%my_rows_idx)) Then
                        If (Any(j_glob .eq. system%overlap%my_cols_idx)) Then
                           i_loc = glob2loc(system%overlap%my_rows_idx,i_glob)
                           j_loc = glob2loc(system%overlap%my_cols_idx,j_glob)
#endif
                           zt = dcmplx (ololo(ilo1, ilo2, ias), 0.0)
                           Call Hermitianmatrix_indexedupdate (system%overlap, j_loc, i_loc, zt)
#ifdef MPI
write (*,*) MPIglobal%rank, 'glob', i_glob, j_glob, 'loc', i_loc, j_loc, 'value', zt
                        End If
                     End If
#endif
                  End If
               End Do
            End If
         End Do
      End Do
      Return
End Subroutine


