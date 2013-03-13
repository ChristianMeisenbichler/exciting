!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
#ifdef MPI      
Subroutine hmlalon (hamilton, is, ia, ngp, apwalm, ngp_loc_row, ngp_loc_col, cols_loc2glob)
#else
Subroutine hmlalon (hamilton, is, ia, ngp, apwalm)
#endif
      Use modinput,        Only: input
      Use mod_muffin_tin,  Only: idxlm, lmmaxapw
      Use mod_APW_LO,      Only: apword, apwordmax, nlorb, lorbl
      Use mod_atoms,       Only: natmtot, idxas
      Use mod_eigensystem, Only: gntyry, idxlo, hloa
      Use modfvsystem,     Only: HermitianMatrix, Hermitianmatrix_indexedupdate
      Implicit None
! arguments
      Type (HermitianMatrix), Intent (Inout) :: hamilton
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
#ifdef MPI      
      Complex (8), Intent (In) :: apwalm (ngp_loc_row, apwordmax, lmmaxapw, natmtot) !SPLIT first dimension over procs
      Integer, Intent (In)     :: ngp_loc_row, ngp_loc_col !TODO better name
      Integer, Dimension(hamilton%ncols_loc), Intent (In) :: cols_loc2glob
#else 
      Complex (8), Intent (In) :: apwalm (ngp, apwordmax, lmmaxapw, natmtot) 
#endif
!
! local variables
      Integer :: ias, io, ilo, i_glob, i_loc, j
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Complex (8) zsum, zt1

      
      ias = idxas (ia, is)
#ifdef MPI
      i_loc = ngp_loc_col
#endif
      Do ilo = 1, nlorb (is)
         l1 = lorbl (ilo, is)
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            i_glob = ngp + idxlo (lm1, ilo, ias)
#ifdef MPI
            if (Any(cols_loc2glob .eq. i_glob)) then
               i_loc = i_loc+1 !!! assuming that the indexing increases
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
#ifdef MPI
                           Do j = 1, ngp_loc_row
                              zt1 = zsum * apwalm (j, io, lm3, ias)
                              hamilton%za(j,i_loc) = hamilton%za(j,i_loc) + conjg(zt1)
                           End Do
#else
                           Do j = 1, ngp
                              zt1 = zsum * apwalm (j, io, lm3, ias)
                              Call Hermitianmatrix_indexedupdate(hamilton, i_glob, j, conjg(zt1))
                           End Do
#endif
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
