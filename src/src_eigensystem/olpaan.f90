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
#ifdef MPI
Subroutine olpaan (overlap, is, ia, ngp, apwalm, ngp_loc)
#else
Subroutine olpaan (overlap, is, ia, ngp, apwalm)
#endif
! !USES:
      Use modinput,        Only: input
      Use mod_muffin_tin,  Only: idxlm, lmmaxapw, lmmaxmat
      Use mod_APW_LO,      Only: apword, apwordmax
      Use mod_atoms,       Only: natmtot, idxas
      Use modfvsystem,     Only: HermitianMatrix, ComplexMatrix, &
                                 newmatrix, deletematrix, &
                                 HermitianMatrixMatrix
#ifdef MPI
     Use modmpi,           Only: DISTRIBUTE_COLS
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!   overlap  : Overlap matrix to update (inout,HermitianMatrix)
!   is       : species number (in,integer)
!   ia       : atom number (in,integer)
!   ngp      : number of G+p-vectors (in,integer)
!   apwalm   : APW matching coefficients
!              (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!              when using MPI it is distributed along the first dimension
!   ngp_loc  : size (1st dimension) of local part of apwalm
!
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the overlap matrix.
!   When MPI is used the overlap matrix is expected to have distributed columns
!
! !REVISION HISTORY:
!   Parallelized and with new matrix types March 2013 (G. Huhs - BSC)
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Type (HermitianMatrix), Intent (Inout) :: overlap
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
#ifdef MPI
      Complex (8), Intent (In) :: apwalm (ngp_loc, apwordmax, lmmaxapw, &
     & natmtot)   !SPLIT first dimension over procs
      Integer, Intent (In) :: ngp_loc
#else
      Complex (8), Intent (In) :: apwalm (ngp, apwordmax, lmmaxapw, &
     & natmtot)   
#endif
!
! local variables
      Integer :: ias, l, m, lm, io,naa
      Type(ComplexMatrix) :: zm1

#ifdef MPI
      Call newmatrix(zm1, lmmaxmat*apwordmax, ngp, DISTRIBUTE_COLS)
#else
      Call newmatrix(zm1, lmmaxmat*apwordmax, ngp)
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

      call HermitianMatrixMatrix(overlap,zm1,zm1,lmmaxmat*apwordmax,naa,ngp)

      Call deletematrix(zm1)

      Return
End Subroutine
