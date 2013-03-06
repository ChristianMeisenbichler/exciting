!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
!
!
Subroutine seceqnfv (nmatp, ngp, igpig, vgpc, apwalm, evalfv, evecfv)
  ! !USES:
      Use modfvsystem,     Only: evsystem, newsystem
      Use mod_muffin_tin,  Only: lmmaxapw
      Use mod_APW_LO,      Only: apwordmax
      Use mod_Gkvector,    Only: ngkmax
      Use mod_atoms,       Only: natmtot
      Use mod_eigenvalue_occupancy, Only: nstfv
      Use mod_eigensystem, Only: nmatmax
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
  !   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
  !   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
  !   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  !   Solves the secular equation,
  !   $$ (H-\epsilon O)b=0, $$
  !   for the all the first-variational states of the input $k$-point.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: nmatp
      Integer, Intent (In) :: ngp
      Integer, Intent (In) :: igpig (ngkmax)
      Real (8), Intent (In) :: vgpc (3, ngkmax)
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Real (8), Intent (Out) :: evalfv (nstfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv)
  ! local variables
      Type (evsystem) :: system

  ! allocatable arrays

!
  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!
!

!distribute the relevant part of apwalm 


      Call newsystem (system,   nmatp)
      Call hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
!
  !------------------------------------!
  !     solve the secular equation     !
  !------------------------------------!
     Call solvewithlapack(system,nstfv,evecfv,evalfv)

End Subroutine seceqnfv
!EOC
