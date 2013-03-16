!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hmllolon (system, is, ia)
      Use modinput,        Only: input
      Use mod_muffin_tin,  Only: idxlm
      Use mod_APW_LO,      Only: lorbl, nlorb
      Use mod_atoms,       Only: idxas
      Use mod_eigensystem, Only: gntyry, hlolo, idxlo
      Use modfvsystem,     Only: evsystem, Hermitianmatrix_indexedupdate
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
      Integer :: ias, ilo1, ilo2
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Complex (8) zsum

#ifdef MPI
      i_loc = system%ngp_loc_rows
#else
      i_loc => i_glob
      j_loc => j_glob
#endif
      ias = idxas(ia, is)
      Do ilo1 = 1, nlorb (is)
         l1 = lorbl (ilo1, is)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            i_glob = system%ngp + idxlo(lm1, ilo1, ias)
#ifdef MPI
            if (Any(i_glob .eq. system%hamilton%my_rows_idx)) then
               i_loc = i_loc+1 !!! assuming that the indexing increases
                                 ! => problem if indexing or order of loops changes
                                 ! an array linking from global to local indices would be better
               j_loc = system%ngp_loc_cols
#endif

               Do ilo2 = 1, nlorb (is)
                  l3 = lorbl (ilo2, is)
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                     j_glob = system%ngp + idxlo(lm3, ilo2, ias)
#ifdef MPI
                     if (Any(j_glob .eq. system%hamilton%my_cols_idx)) then
                        j_loc = j_loc+1 !!! assuming that the indexing increases
                                          ! => problem if indexing or order of loops changes
                                          ! an array linking from global to local indices would be better
#endif
                        If (i_glob .Le. j_glob) Then
                           zsum = 0.d0
                           Do l2 = 0, input%groundstate%lmaxvr
                              If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                                 Do m2 = - l2, l2
                                    lm2 = idxlm (l2, m2)
                                    zsum = zsum + gntyry (lm1, lm2, lm3) * &
                                 & hlolo (ilo1, ilo2, lm2, ias)
                                 End Do
                              End If
                           End Do

! calculate the matrix elements
                           Call Hermitianmatrix_indexedupdate (system%hamilton, j_loc, i_loc, zsum)
                        End If
#ifdef MPI
                     End If
#endif
                  End Do
               End Do
#ifdef MPI
            End If
#endif
         End Do
      End Do
      Return
End Subroutine
