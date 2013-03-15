!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine hmlalon (system, is, ia, apwalm)
      Use modinput,        Only: input
      Use mod_muffin_tin,  Only: idxlm, lmmaxapw
      Use mod_APW_LO,      Only: apword, apwordmax, nlorb, lorbl
      Use mod_atoms,       Only: natmtot, idxas
      Use mod_eigensystem, Only: gntyry, idxlo, hloa
      Use modfvsystem,     Only: evsystem, Hermitianmatrix_indexedupdate
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Complex (8), Intent (In) :: apwalm (system%ngp_loc_rows, apwordmax, lmmaxapw, natmtot) !SPLIT first dimension over procs
!
! local variables
#ifdef MPI
      Integer :: j_loc
      Integer :: j_glob
#else
      Integer, Pointer :: j_loc
      Integer, Target  :: j_glob
#endif
      Integer :: ias, io, ilo, i
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Complex (8) zsum, zt1

      
      ias = idxas (ia, is)
#ifdef MPI
      j_loc = system%ngp_loc_cols
#else
      j_loc => j_glob
#endif
      Do ilo = 1, nlorb (is)
         l1 = lorbl (ilo, is)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            j_glob = system%ngp + idxlo (lm1, ilo, ias)
#ifdef MPI
            if (Any(j_glob .eq. system%hamilton%my_cols_idx)) then
               j_loc = j_loc+1 !!! assuming that the indexing increases
                               ! => problem if indexing or order of loops changes
                               ! an array linking between global and local indices would be better
#endif
               Do l3 = 0, input%groundstate%lmaxmat
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                     Do io = 1, apword (l3, is)

                        zsum = 0.d0
                        Do l2 = 0, input%groundstate%lmaxvr
                           If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                              Do m2 = - l2, l2
                                 lm2 = idxlm (l2, m2)
                                 zsum = zsum + gntyry (lm1, lm2, lm3) * &
                                 & hloa (ilo, io, l3, lm2, ias)
                              End Do
                           End If
                        End Do

                        ! calculate the matrix elements
                        If (Abs(dble(zsum))+Abs(aimag(zsum)) .Gt. 1.d-20) Then
                           Do i = 1, system%ngp_loc_rows
                              zt1 = zsum * apwalm (i, io, lm3, ias)
                              Call Hermitianmatrix_indexedupdate(system%hamilton, j_loc, i, conjg(zt1))
!commented: performance-version for MPI
! #ifdef MPI
!                               system%hamilton%za(i,j_loc) = system%hamilton%za(i,j_loc) + conjg(zt1)
! #else
!                               Call Hermitianmatrix_indexedupdate(system%hamilton, j_glob, i, conjg(zt1))
! #endif
                           End Do
                        End If
                     End Do
                  End Do
               End Do
#ifdef MPI
            End If
#endif
         End Do
      End Do
      Return
End Subroutine
