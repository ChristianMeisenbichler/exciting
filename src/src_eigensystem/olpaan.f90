!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: olpaan
! !INTERFACE:
!
!
Subroutine olpaan (system, is, ia, apwalm)
! !USES:
      Use modinput,        Only: input
      Use mod_muffin_tin,  Only: idxlm, lmmaxapw, lmmaxmat
      Use mod_APW_LO,      Only: apword, apwordmax
      Use mod_atoms,       Only: natmtot, idxas
      Use modfvsystem,     Only: evsystem, ComplexMatrix, &
                                 newmatrix, deletematrix, &
                                 HermitianMatrixMatrix
#ifdef MPI
     Use modmpi,           Only: DISTRIBUTE_COLS
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!   system   : EVSystem with the Overlap matrix to update (inout,evsystem)
!   is       : species number (in,integer)
!   ia       : atom number (in,integer)
!   apwalm   : APW matching coefficients
!              (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!              when using MPI it is distributed along the first dimension
!
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the overlap matrix.
!   When MPI is used the overlap matrix is expected to have distributed columns
!   and apwalm's first dimension has to be distributed according to that
!
! !REVISION HISTORY:
!   Parallelized and with new matrix types March 2013 (G. Huhs - BSC)
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Type (evsystem), Intent (Inout) :: system
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Complex (8), Intent (In) :: apwalm (system%ngp_loc_cols, apwordmax, lmmaxapw, &
     & natmtot)   !SPLIT first dimension over procs
!
! local variables
      Integer :: ias, l, m, lm, io,naa
      Type(ComplexMatrix) :: zm1

#ifdef MPI
      Call newmatrix(zm1, lmmaxmat*apwordmax, system%ngp, DISTRIBUTE_COLS)
#else
      Call newmatrix(zm1, lmmaxmat*apwordmax, system%ngp)
#endif

      ias = idxas (ia, is)
      naa=0

      Do l = 0, input%groundstate%lmaxmat
         Do m = - l, l
            lm = idxlm (l, m)
            Do io = 1, apword (l, is)
               naa=naa+1
               zm1%za(naa,:)=apwalm(1:zm1%ncols_loc, io, lm, ias)
            End Do
         End Do
      End Do

      call HermitianMatrixMatrix(system%overlap,zm1,zm1,lmmaxmat*apwordmax,naa,system%ngp)

      Call deletematrix(zm1)

      Return
End Subroutine
