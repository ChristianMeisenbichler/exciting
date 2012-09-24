!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqn
!
!
Subroutine seceqn (ik, evalfv, evecfv, evecsv,flags)
  ! !USES:
      Use modinput
      Use modmain
      Use modmpi
      Use sclcontroll
      Use diisinterfaces
      use mod_libapw
!
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  !   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
  !   flags  : flag bits used by lapw_exec. 1:evec,2evec+density
  ! !DESCRIPTION:
  !   Solves the first- and second-variational secular equations. See routines
  !   {\tt match}, {\tt seceqnfv}, {\tt seceqnss} and {\tt seceqnsv}.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik,flags
      Real (8), Intent (Out) :: evalfv (nstfv, nspnfv)
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
  ! local variables
      Integer :: ispn
!
!
  ! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)
#ifdef KSMP
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO
#endif
  !
  !-IMPORTANT: the first-variational spinor index and the k-point index have been
  ! swapped in the following arrays: ngk, igkig, vgkl, vgkc, gkc, tpgkc, sfacgk
  !
  
#ifdef _LIBAPW_
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
   &sfacgk(:,:,1,ik),apwalm(:,:,:,:,1))
call lapw_execute(ik,apwalm,evalsv(1,ik),occsv(1,ik),densmt,densir,flags)
deallocate(apwalm)
return
#endif

      Do ispn = 1, nspnfv
     ! find the matching coefficients
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
     ! solve the first-variational secular equation
         If (doDIIScycle()) Then
            Call DIISseceqnfv (ik, ispn, apwalm(:, :, :, :, ispn), &
           & vgkc(:, :, ispn, ik), evalfv, evecfv)
!
            If (ik .Eq. lastk(rank)) diiscounter = diiscounter + 1
!
         Else If (doARPACKiteration()) Then
            Call iterativearpacksecequn (ik, ispn, apwalm(1, 1, 1, 1, &
           & ispn), vgkc(1, 1, ispn, ik), evalfv, evecfv)
         Else If (doLAPACKsolver()) Then
            If (tseqit) Then
           ! iteratively
               Call seceqnit (nmat(ispn, ik), ngk(ispn, ik), igkig(:, &
              & ispn, ik), vkl(:, ik), vgkl(:, :, ispn, ik), vgkc(:, :, &
              & ispn, ik), apwalm(:, :, :, :, ispn), evalfv(:, ispn), &
              & evecfv(:, :, ispn))
            Else
           ! directly
               Call seceqnfv (nmat(ispn, ik), ngk(ispn, ik), igkig(:, &
              & ispn, ik), vgkc(:, :, ispn, ik), apwalm(:, :, :, :, &
              & ispn), evalfv(:, ispn), evecfv(:, :, ispn))
            End If
         Else If (.True.) Then
            Write (*,*) "error in solverselect secequn.F90"
!
         End If
      End Do
#ifdef KSMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
      If (isspinspiral()) Then
     ! solve the spin-spiral second-variational secular equation
         Call seceqnss (ik, apwalm, evalfv, evecfv, evecsv)
      Else
     ! solve the second-variational secular equation
         Call seceqnsv (ik, apwalm, evalfv, evecfv, evecsv)
      End If
!
      Deallocate (apwalm)
      Return
End Subroutine seceqn
!EOC
