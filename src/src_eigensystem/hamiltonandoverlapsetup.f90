
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
!
!
Subroutine hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
! !USES:
      Use modfvsystem, Only: evsystem
      Use mod_atoms,   Only: nspecies, natoms, natmtot
      Use mod_timing,  Only: time_hmlaan, time_hmlalon, time_hmllolon, & 
                             time_olpaan, time_olpalon, time_olplolon, &
                             time_hmlistln, time_olpistln, timemat
      Use mod_muffin_tin, Only: lmmaxapw  !only needed for shape declaration of apwalm
      Use mod_APW_LO,     Only: apwordmax !only needed for shape declaration of apwalm
      Use mod_gkvector,   Only: ngkmax    !only needed for shape declaration of apwalm
!
! !INPUT/OUTPUT PARAMETERS:
!   system : 
!   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
! !DESCRIPTION:
!   Sets up the Hamiltonian and overlap matrices
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Type (evsystem)          :: system
      Integer, Intent (In)     :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Integer, Intent (In)     :: igpig (ngkmax)
      Real (8), Intent (In)    :: vgpc (3, ngkmax)
! local variables
      Real (8) :: ts0, ts1
      Integer  :: is, ia
      Real (8) :: cpu0, cpu1
      Real (8) :: threshold
!----------------------------------------!
!     Hamiltonian and overlap set up     !
!----------------------------------------!
!
!
      Call timesec (cpu0)
! set the matrices to zero
!
! muffin-tin contributions
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
             Call timesec (ts0)
            Call hmlaan (system, is, ia, apwalm)
             Call timesec (ts1)
             time_hmlaan=ts1-ts0+time_hmlaan
             Call timesec (ts0)
            Call hmlalon (system, is, ia, apwalm)
             Call timesec (ts1)
             time_hmlalon=ts1-ts0+time_hmlalon
             Call timesec (ts0)
            Call hmllolon (system%hamilton, is, ia, ngp)
             Call timesec (ts1)
             time_hmllolon=ts1-ts0+time_hmllolon
             Call timesec (ts0)
            Call olpaan (system, is, ia, apwalm)
             Call timesec (ts1)
             time_olpaan=ts1-ts0+time_olpaan
             Call timesec (ts0)
            Call olpalon (system, is, ia, apwalm)
             Call timesec (ts1)
             time_olpalon=ts1-ts0+time_olpalon
             Call timesec (ts0)
            Call olplolon (system%overlap, is, ia, ngp)
             Call timesec (ts1)
             time_olplolon=ts1-ts0+time_olplolon
         End Do
      End Do
!
! interstitial contributions
       Call timesec (ts0)
      Call hmlistln (system, igpig, vgpc)
       Call timesec (ts1)
       time_hmlistln=ts1-ts0+time_hmlistln
       Call timesec (ts0)
      Call olpistln (system%overlap, ngp, igpig)
       Call timesec (ts1)
       time_olpistln=ts1-ts0+time_olpistln
      threshold = 1e-16
!call HermitianMatrixTruncate(system%hamilton,threshold)
!call HermitianMatrixTruncate(system%overlap,threshold)
!
!
!

#ifdef DEBUGHO
      Write (*,*) "apwalm", apwalm
      prefix = "H"
      Call HermitianMatrixToFiles (system%hamilton, prefix)
      prefix = "O"
      Call HermitianMatrixToFiles (system%overlap, prefix)
      Write (*,*) "wrote"
      Stop
#endif
!
      Call timesec (cpu1)
      timemat = timemat + cpu1 - cpu0

End Subroutine
